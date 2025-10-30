rm(list = ls())
library(dplyr)
library(tidyr)
library(Matrix)

load("./data/california_smoking.rda")
california_smoking

Rcpp::sourceCpp("./R/20_mcmc.cpp")
Rcpp::sourceCpp("./R/40_geweke_latest.cpp")
source("R/01_utils.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/10_sc_spillover.R")
source("R/03_utils_plot.R")
source("R/40_geweke_latest.R")
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


##########################################################################################
##########################################################################################
# Prior Sensitivity
##########################################################################################
##########################################################################################
## ---- 事後MCMC（観測データ固定） ----
run_mcmc_for_posterior <- function(
  Yc_obs,               # T0 x N（観測データ: 前期コントロール行列）
  W_raw, w_raw,
  alpha_hat_scaled,     # 長さN（BSCMのαを列標準化に整合）
  Xc_pre = NULL,        # T0 x N x K あるいは NULL
  p = 0L,               # 因子数（まずは 0 を推奨）
  a0 = 1.0, b0 = 1.0,   # σ² ~ IG(a0,b0)
  rho_support = NULL,   # c(lo,hi); NULLなら固有値境界から計算
  step_rho = 0.05,
  M_burn = 5000L, M_keep = 20000L,
  thin = 1L,            # 事後サンプルの間引き（1=無間引き）
  seed = 123
){
  stopifnot(is.matrix(Yc_obs))
  set.seed(seed)

  T0 <- nrow(Yc_obs); N <- ncol(Yc_obs)
  K  <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  # ---- X の三次元化（NULL→空配列）----
  X_use <- if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else {
    stopifnot(dim(Xc_pre)[1] == T0, dim(Xc_pre)[2] == N)
    Xc_pre
  }

  # ---- (w | W) の正規化 ----
  nm <- normalize_joint_wW(W_raw, w_raw)
  W  <- nm$W_new; w <- nm$w_new

  # ---- rho の支持域 ----
  if (is.null(rho_support)) {
    bnd <- compute_bnd(W, c_stability = 0.95)
    rho_lo <- -bnd; rho_hi <-  bnd
  } else {
    stopifnot(length(rho_support) == 2L, rho_support[1] < rho_support[2])
    rho_lo <- as.numeric(rho_support[1]); rho_hi <- as.numeric(rho_support[2])
  }

  # ---- 初期状態（事前整合）----
  rinvgamma1 <- function(a, b) 1/rgamma(1, shape=a, rate=b)
  state <- list(
    rho    = as.numeric(runif(1, rho_lo, rho_hi)),
    sigma2 = as.numeric(rinvgamma1(a0, b0)),
    beta   = if (K>0) rnorm(K, 0, 1) else numeric(0),
    Eta    = if (p>0) matrix(0, N, p) else matrix(0, N, 0),
    Gamma  = if (p>0) matrix(0, p, T0) else matrix(0, 0, T0)
  )

  # ---- burn-in：観測固定で 1 ステップカーネルのみ ----
  for (m in seq_len(M_burn)) {
    state <- scspill_one_step_cpp(
      Yc_data = Yc_obs,
      W_use = W, w_use = w,
      alpha_hat_scaled = alpha_hat_scaled,
      T0 = T0, N = N, Xc_pre = X_use, K = K, p = p,
      state_in = state, a0 = a0, b0 = b0, step_rho = step_rho,
      rho_lo = rho_lo, rho_hi = rho_hi
    )
  }

  # ---- keep：間引き対応 ----
  keep_indices <- seq_len(M_keep * thin)
  draws <- vector("list", length = M_keep)
  k <- 0L
  for (m in keep_indices) {
    state <- scspill_one_step_cpp(
      Yc_data = Yc_obs,
      W_use = W, w_use = w,
      alpha_hat_scaled = alpha_hat_scaled,
      T0 = T0, N = N, Xc_pre = X_use, K = K, p = p,
      state_in = state, a0 = a0, b0 = b0, step_rho = step_rho,
      rho_lo = rho_lo, rho_hi = rho_hi
    )
    if (m %% thin == 0L) {
      k <- k + 1L
      draws[[k]] <- state
    }
  }

  # ---- 事後要約 ----
  rho_vec <- vapply(draws, function(s) s$rho,    numeric(1))
  s2_vec  <- vapply(draws, function(s) s$sigma2, numeric(1))

  out_df <- data.frame(
    param = c("rho","sigma2"),
    mean  = c(mean(rho_vec), mean(s2_vec)),
    sd    = c(sd(rho_vec),   sd(s2_vec)),
    q025  = c(unname(quantile(rho_vec, 0.025)), unname(quantile(s2_vec, 0.025))),
    q975  = c(unname(quantile(rho_vec, 0.975)), unname(quantile(s2_vec, 0.975)))
  )

  if (K > 0) {
    # beta はベクトル長 K。列方向に要約
    B <- do.call(cbind, lapply(draws, `[[`, "beta")) # K x M_keep
    beta_summ <- data.frame(
      param = paste0("beta[", seq_len(K), "]"),
      mean  = rowMeans(B),
      sd    = apply(B, 1, sd),
      q025  = apply(B, 1, function(x) quantile(x, 0.025)),
      q975  = apply(B, 1, function(x) quantile(x, 0.975))
    )
    out_df <- rbind(out_df, beta_summ)
  }

  list(theta = out_df,
       raw = list(rho = rho_vec, sigma2 = s2_vec),
       meta = list(M_burn = M_burn, M_keep = M_keep, thin = thin,
                   rho_lo = rho_lo, rho_hi = rho_hi, step_rho = step_rho))
}

## ---- 事前感度ラッパー ----
prior_sensitivity <- function(
  Yc_obs, W_raw, w_raw,
  alpha_hat_scaled,
  Xc_pre = NULL, p = 0L,
  grid,                       # data.frame: a0, b0, rho_lo, rho_hi, step_rho
  M_burn = 5000L, M_keep = 20000L,
  thin = 1L
){
  stopifnot(is.data.frame(grid))
  req <- c("a0","b0","rho_lo","rho_hi","step_rho")
  stopifnot(all(req %in% names(grid)))

  out_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    a0 <- grid$a0[i]; b0 <- grid$b0[i]
    rho_lo <- grid$rho_lo[i]; rho_hi <- grid$rho_hi[i]
    step_rho <- grid$step_rho[i]

    res <- run_mcmc_for_posterior(
      Yc_obs = Yc_obs,
      W_raw = W_raw, w_raw = w_raw,
      alpha_hat_scaled = alpha_hat_scaled,
      Xc_pre = Xc_pre, p = p,
      a0 = a0, b0 = b0,
      rho_support = c(rho_lo, rho_hi),
      step_rho = step_rho,
      M_burn = M_burn, M_keep = M_keep,
      thin = thin,
      seed = 1000 + i
    )
    out_list[[i]] <- list(set = grid[i, , drop = FALSE],
                          theta = res$theta)
  }

  theta_table <- do.call(
    rbind,
    lapply(out_list, function(one) cbind(one$set, one$theta))
  )
  rownames(theta_table) <- NULL

  list(results = out_list, theta_table = theta_table)
}

grid <- data.frame(
  a0 = c(3,5,2),
  b0 = c(1,1,0.5),
  rho_lo = c(-0.5,-0.3,-0.7),
  rho_hi = c( 0.5, 0.3, 0.7),
  step_rho = c(0.05, 0.03, 0.07)
)

Yc_pre_obs <- panel_df %>%
  filter(
    (state != "California") & (year < 1988)
  ) %>%
  pivot_wider(
    id_cols = year,
    names_from = state,
    values_from = cigsale
  ) %>%
  select(-year) %>%
  as.matrix()

Y0_pre <- panel_df %>%
  filter(
    (state == "California") & (year < 1988)
  ) %>%
  select(cigsale) %>%
  as.matrix()

Xc_pre <- panel_df %>%
  filter(
    (state != "California") & (year < 1988)
  ) %>%
  pivot_wider(
    id_cols = year,
    names_from = state,
    values_from = retprice
  ) %>%
  select(-year) %>%
  as.matrix()

dim(Xc_pre) <- c(nrow(Xc_pre), ncol(Xc_pre), 1)

alpha_hat <- colMeans(fit$alpha_draws)
sds_Yc <- apply(Yc_pre_obs, 2, sd)
sds_Yc[!is.finite(sds_Yc) | sds_Yc < 1e-8] <- 1
alpha_hat_scaled <- alpha_hat / sds_Yc

sens <- prior_sensitivity(
  Yc_obs = Yc_pre_obs,       # 観測（T0 x N）
  W_raw = W, w_raw = w,
  alpha_hat_scaled = alpha_hat_scaled,
  Xc_pre = Xc_pre, p = 0L,
  grid = grid,
  M_burn = 50000L, M_keep = 100000L
)

sens$theta_table   # 設定別の事後要約



##########################################################################################
# Prior Predictive
##########################################################################################
## ---- Rだけでの前進シミュレータ（β と 因子を反映） ----
simulate_Yc_forward_R <- function(T0, W_use, w_use, alpha_hat_scaled,
                                  rho, sigma2, Xc_pre, beta, Eta, Gamma) {
  N <- nrow(W_use); K <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]
  A <- W_use + w_use * as.numeric(alpha_hat_scaled)
  A <- matrix(A, nrow=N, ncol=N)  # 明示
  I <- diag(N)
  Yc <- matrix(NA_real_, T0, N)
  sd_eps <- sqrt(max(sigma2, 1e-12))

  for (t in seq_len(T0)) {
    mu <- rep(0, N)
    if (K > 0) {
      arr <- Xc_pre[t, , , drop = FALSE]
      Nloc <- dim(arr)[2]; Kloc <- dim(arr)[3]
      Xt <- matrix(arr, nrow = Nloc, ncol = Kloc)
      mu <- mu + as.numeric(Xt %*% beta)
    }
    if (ncol(Eta) > 0) mu <- mu + as.numeric(Eta %*% Gamma[,t])
    eps <- rnorm(N, 0, sd_eps)
    rhs <- mu + eps
    # (I - rho A) Y_t = rhs
    M <- I - rho * A
    # 数値安定のための fallback
    Yt <- tryCatch(solve(M, rhs), error=function(e) {
      AtA <- crossprod(M)
      AtA[cbind(seq_len(N), seq_len(N))] <- AtA[cbind(seq_len(N), seq_len(N))] + 1e-10
      solve(AtA, crossprod(M, rhs))
    })
    Yc[t,] <- Yt
  }
  Yc
}

## ---- 事前分布からの完全サンプラ ----
prior_sampler_theta_full <- function(T0, N, K, p, a0, b0,
                                     W_use, w_use, alpha_hat_scaled,
                                     rho_support = NULL) {
  if (is.null(rho_support)) {
    bnd <- compute_bnd(W_use, c_stability = 0.95)
    rho <- runif(1, -bnd, bnd)
  } else {
    rho <- runif(1, rho_support[1], rho_support[2])
  }
  sigma2 <- 1/rgamma(1, shape=a0, rate=b0)
  beta   <- if (K>0) rnorm(K, 0, 1) else numeric(0)
  Eta    <- if (p>0) matrix(rnorm(N*p, 0, 1), N, p) else matrix(0, N, 0)
  Gamma  <- if (p>0) matrix(rnorm(p*T0, 0, 1), p, T0) else matrix(0, 0, T0)
  list(rho=rho, sigma2=sigma2, beta=beta, Eta=Eta, Gamma=Gamma)
}

## ---- 比較統計（必要に応じて追加）----
ppc_stats <- function(Yc, Y0_pre, W_use, w_use){
  yc <- as.numeric(Yc)
  N <- ncol(Yc); T0 <- nrow(Yc)
  wyc <- as.numeric(Yc %*% as.numeric(w_use))
  spatial_q <- sum(diag(Yc %*% W_use %*% t(Yc))) / (N*T0)
  ac1 <- tryCatch({
    Yd <- scale(t(Yc), center=TRUE, scale=FALSE)  # N x T0
    num <- rowSums(Yd[, -1, drop=FALSE] * Yd[, -ncol(Yd), drop=FALSE])
    den <- rowSums(Yd * Yd)
    mean(num/den, na.rm=TRUE)
  }, error=function(e) NA_real_)
  c(
    yc_mean = mean(yc),
    log_yc_var = log(var(yc) + 1e-12),
    spatial_quadratic = spatial_q,
    corr_y0_wyc = if (sd(Y0_pre)>0 && sd(wyc)>0) cor(Y0_pre, wyc) else NA_real_,
    ac1 = ac1
  )
}

## ---- Prior Predictive 本体 ----
prior_predictive <- function(
  Y0_pre, Yc_obs = NULL,        # 観測統計線を引く場合は Yc_obs を与える
  W_raw, w_raw, alpha_hat_scaled,
  Xc_pre = NULL, p = 0L,
  a0=3, b0=1, rho_support=NULL,
  R=2000L, seed=123
){
  set.seed(seed)
  T0 <- length(Y0_pre)
  N  <- nrow(W_raw)
  K  <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  nm <- normalize_joint_wW(W_raw, w_raw)
  W  <- nm$W_new; w <- nm$w_new

  # 観測統計（任意）
  obs_stat <- if (!is.null(Yc_obs)) ppc_stats(Yc_obs, Y0_pre, W, w) else NULL

  stat_mat <- matrix(NA_real_, R, 5)
  colnames(stat_mat) <- c("yc_mean","log_yc_var","spatial_quadratic","corr_y0_wyc","ac1")

  for (r in seq_len(R)) {
    th <- prior_sampler_theta_full(
      T0, N, K, p, a0, b0, W, w, alpha_hat_scaled, rho_support
    )
    Yc_sim <- simulate_Yc_forward_R(
      T0, W, w, alpha_hat_scaled,
      th$rho, th$sigma2,
      if (is.null(Xc_pre)) array(0, c(T0,N,0)) else Xc_pre,
      th$beta, th$Eta, th$Gamma
    )
    stat_mat[r, ] <- ppc_stats(Yc_sim, Y0_pre, W, w)
  }

  list(stat = as.data.frame(stat_mat),
       observed = obs_stat)
}

# 観測統計と prior predictive の分布を比較
ppc <- prior_predictive(
  Y0_pre = Y0_pre,
  Yc_obs = Yc_pre_obs,
  W_raw = W, w_raw = w,
  alpha_hat_scaled = alpha_hat_scaled,
  Xc_pre = Xc_pre, p = 0L,
  a0 = 5, b0 = 1, rho_support = c(-0.3, 0.3),
  R = 20000
)

# 例: 分布に観測値を重ねて可視化
op <- par(mfrow=c(2,3))
for (nm in colnames(ppc$stat)) {
  hist(ppc$stat[[nm]], breaks=30, main=paste("Prior predictive:", nm),
       xlab=nm, col="grey90", border="white")
  if (!is.null(ppc$observed) && !is.na(ppc$observed[nm])) {
    abline(v = ppc$observed[nm], col=2, lwd=2)
  }
}
par(op)













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


# --- dimensions ---
T0 <- 15
N  <- 8
K  <- 2     # Xを入れる
p  <- 1     # 因子も入れる

# --- W, w, alpha_hat_scaled の用意 ---
W_raw <- matrix(0, N, N)
for (i in 1:N) {
  if (i > 1) W_raw[i, i-1] <- 1
  if (i < N) W_raw[i, i+1] <- 1
}
w_raw <- rep(0, N); w_raw[1] <- 1

normalize_joint_wW <- function(W, w){
  stopifnot(nrow(W)==length(w), ncol(W)==nrow(W))
  joint <- cbind(w, W)
  rs <- rowSums(joint); rs[rs==0] <- 1
  joint_norm <- joint/rs
  list(W_new = joint_norm[,-1,drop=FALSE],
       w_new = joint_norm[,1,drop=TRUE])
}
normed <- normalize_joint_wW(W_raw, w_raw)
W_use <- normed$W_new; w_use <- normed$w_new

# α の基準スケール：ここでは単純に固定値（本番はBSCMのスケールと整合させる）
alpha_hat_scaled <- rnorm(N, 0, 0.4)

# --- 共変量 X (T0 x N x K) の生成 ---
Xc_pre <- array(0, dim=c(T0, N, K))
for (t in 1:T0) for (i in 1:N) for (k in 1:K)
  Xc_pre[t,i,k] <- rnorm(1)

# --- 真のパラメータ（検証用） ---
rho_true    <- 0.25
sigma2_true <- 0.5^2
beta_true   <- c(0.7, -0.5)
Eta_true    <- matrix(rnorm(N*p, 0, 0.5), N, p)
Gamma_true  <- matrix(0, p, T0)
phi_true    <- 0.5
s2g_true    <- 0.2
for (t in 1:T0){
  if (t==1) Gamma_true[,t] <- rnorm(p, 0, sqrt(s2g_true/(1-phi_true^2)))
  if (t>=2) Gamma_true[,t] <- phi_true*Gamma_true[,t-1] + rnorm(p, 0, sqrt(s2g_true))
}

# --- 事前（JDT用） ---
a0 <- 3.0; b0 <- 1.0    # IG(a0,b0) for sigma^2
step_rho <- 0.05

# --- Rcpp 側の関数をロード（事前に sourceCpp しておく）
# Rcpp::sourceCpp("40_geweke_full.cpp")   # ここに simulate_Yc_forward_cpp/scspill_one_step_cpp がある想定

# --- MC側での forward 生成と同型のシミュレータをRで（チェック用） ---
sim_forward_R <- function(T0, W_use, w_use, alpha_hat_scaled,
                          rho, sigma2, Xc_pre, beta, Eta, Gamma){
  N <- nrow(W_use); K <- dim(Xc_pre)[3]
  A <- W_use + w_use %*% t(alpha_hat_scaled)
  I <- diag(N)
  Yc <- matrix(NA_real_, T0, N)
  for (t in 1:T0){
    mu <- rep(0, N)
    if (K>0) {
      Xt <- Xc_pre[t,,]     # N x K
      mu <- mu + as.vector(Xt %*% beta)
    }
    if (ncol(Eta)>0) mu <- mu + as.vector(Eta %*% Gamma[,t])
    eps <- rnorm(N, 0, sqrt(sigma2))
    rhs <- mu + eps
    M <- I - rho * A
    Yc[t,] <- as.vector(solve(M, rhs))
  }
  Yc
}

# --- テストで使う観測 treated pre-outcome（相関統計用）---
Y0_pre <- rnorm(T0)

# --- JDT ラッパ（既存の geweke_jdt_full を呼ぶ想定） ---
#   ※ 前進シミュレーションで beta/Eta/Gamma を使うこと（C++側修正後）を前提。
#   ※ M1, M2 は小さめでも動作確認はできるが、安定性のため十分大きく。
M1 <- 1000000L
M2 <- 1000000L
burn_in <- 2000000L

# 初期状態ドロー関数は既存のものを使用（beta の初期長さ K を保証）
draw_beta0 <- function(K, s2=1) if (K>0) rnorm(K, 0, sqrt(s2)) else numeric(0)

# 実行（ご自身の geweke_jdt_full を呼んでください。引数名は実装に合わせる）
out_ok <- geweke_jdt_full(
  Y0_pre           = Y0_pre,
  Yc_pre_like_dims = c(T0, N),
  W_raw            = W_use,
  w_raw            = w_use,
  alpha_hat_scaled = alpha_hat_scaled,
  Xc_pre           = Xc_pre,
  p                = p,
  M1               = M1,
  M2               = M2,
  burn_in          = burn_in,
  a0               = a0,
  b0               = b0,
  step_rho         = step_rho,
  g_fn             = default_g_fn,
  verbose          = TRUE
)

print(out_ok$summary)


alpha_hat_scaled_bad <- alpha_hat_scaled * 1.05

out_bad <- geweke_jdt_full(
  Y0_pre           = Y0_pre,
  Yc_dims = c(T0, N),
  W_raw            = W_use,
  w_raw            = w_use,
  alpha_hat_scaled = alpha_hat_scaled_bad,  # ★ SC側ミスマッチ
  Xc_pre           = Xc_pre,
  p                = p,
  M1               = M1,
  M2               = M2,
  burn_in          = burn_in,
  a0               = a0,
  b0               = b0,
  step_rho         = step_rho,
  g_fn             = default_g_fn,
  verbose          = TRUE
)

print(out_bad$summary)
