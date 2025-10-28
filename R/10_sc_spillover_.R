# =========================================================
# sc_spillover(): 二段階推定版
#   Step1: hs_alpha_gibbs_cpp で alpha を推定
#   Step2: sar_full_sampler_cpp_step2 で alpha_hat 固定の下で rho 等を推定
# 効果は alpha_hat を固定し、rho の事後のみ反映
# =========================================================
sc_spillover <- function(
  data,
  treated_unit,
  T0 = NULL,
  w,
  W,
  X = NULL,
  p_factors = 1,
  M = 2000,
  burn = 1000,
  seed = 123,
  verbose = TRUE,
  y = "y",
  unit_col = "unit",
  time_col = "time",
  treatment_dummy,
  step_rho = 0.05,
  step_alpha = 0.05 # 互換性のため残置（本関数では未使用）
) {
  stopifnot(is.data.frame(data))
  if (!is.character(y)) {
    y <- as.character(substitute(y))
  }

  # 必須列
  required_cols <- c(unit_col, time_col, y, treatment_dummy)
  if (!all(required_cols %in% names(data))) {
    stop(sprintf(
      "`data` に列 %s が必要です。",
      paste0(required_cols, collapse = ", ")
    ))
  }

  # Rcpp 実装のロード（未ロードなら）
  if (!exists("hs_alpha_gibbs_cpp") || !exists("sar_full_sampler_cpp_step2")) {
    if (file.exists("20_mcmc.cpp")) {
      Rcpp::sourceCpp("20_mcmc.cpp")
    } else if (file.exists("R/20_mcmc.cpp")) {
      Rcpp::sourceCpp("R/20_mcmc.cpp")
    } else {
      stop(
        "20_mcmc.cpp が見つかりません。`hs_alpha_gibbs_cpp` と `sar_full_sampler_cpp_step2` を含めてください。"
      )
    }
  }

  set.seed(seed)

  # 介入開始の同定 → T0
  times <- sort(unique(data[[time_col]]))
  treat_series <- data[
    data[[unit_col]] == treated_unit,
    c(time_col, treatment_dummy)
  ]
  t_start <- min(treat_series[[time_col]][treat_series[[treatment_dummy]] == 1])
  T0 <- if (is.null(T0)) sum(times < t_start) else as.integer(T0)

  # 前処理（Y0/Yc と X の整形）
  prep <- scspill_prep_X(
    data,
    treated_unit = treated_unit,
    T0 = T0,
    X = X,
    y_col = y,
    unit_col = unit_col,
    time_col = time_col
  )

  # 行標準化（ゼロ割りロバスト実装を想定）
  W <- row_normalize(W)
  # w は Step2 側（C++）で **L2 正規化**されます。ここでは触れません。

  Y0_pre <- prep$Y0_pre
  Yc_pre <- prep$Yc_pre
  Y0_post <- prep$Y0_post
  Yc_post <- prep$Yc_post
  Xc_pre <- prep$Xc_pre
  N <- ncol(Yc_pre)
  T1 <- nrow(Yc_post)

  #--------------------------------------
  # Step 1: BSCM による alpha の推定
  #   Y0_t = alpha' Yc_t + v_t, HS 事前
  #--------------------------------------
  if (verbose) {
    message("[Step 1] Sampling alpha via BSCM (horseshoe prior)...")
  }
  alpha_draws <- hs_alpha_gibbs_cpp(
    Y0_pre = Y0_pre,
    control_outcome_pre = Yc_pre,
    iteration = M,
    burn = burn,
    verbose = verbose
  )
  colnames(alpha_draws) <- colnames(Yc_pre)
  alpha_hat <- colMeans(alpha_draws)

  #--------------------------------------
  # Step 2: \u03B1̂ 固定で SCSPILL の \u03C1 等を推定
  #  (I - rho W - rho w alpha^T) Yc_t = X_t beta + Lambda F_t + u_t
  #--------------------------------------
  if (verbose) {
    message("[Step 2] Sampling rho (and others) with fixed alpha_hat...")
  }
  # X の受け渡し：NULL か、長さ T0*N*K のベクトル
  K <- 0L
  Xvec <- NULL
  if (!is.null(Xc_pre)) {
    if (is.array(Xc_pre) && length(dim(Xc_pre)) == 3L) {
      K <- dim(Xc_pre)[3]
      stopifnot(dim(Xc_pre)[1] == T0, dim(Xc_pre)[2] == N)
      Xvec <- as.numeric(aperm(Xc_pre, c(1, 2, 3)))
    } else {
      # すでにベクトル化されているとみなす
      # K は外側から与えられていないので 0/非0 の整合は scspill_prep_X 側に依存
      Xvec <- as.numeric(Xc_pre)
      K <- length(Xvec) / (T0 * N)
      if (abs(K - round(K)) > 1e-8) {
        stop("Xc_pre の次元が (T0*N*K) に整合しません。")
      }
      K <- as.integer(round(K))
    }
  }

  w_l2 <- as.numeric(w)
  s <- sqrt(sum(w_l2^2))
  if (s > 0) {
    w_l2 <- w_l2 / s
  }

  sar <- sar_full_sampler_cpp_step2(
    Yc_pre = Yc_pre,
    alpha_hat = alpha_hat, # ★ 固定して渡す（C++ 側で Yc スケールと整合）
    Xc_pre_ = if (!is.null(Xvec)) Xvec else R_NilValue,
    T0 = T0,
    N = N,
    K = K,
    p = as.integer(p_factors),
    w_in = as.numeric(w_l2),
    W = as.matrix(W),
    iteration = M,
    burn = burn,
    step_rho = step_rho,
    a0 = 1.0,
    b0 = 1.0,
    verbose = verbose
  )

  rho_draws <- as.numeric(sar$rho)
  rho_hat <- mean(rho_draws)

  #--------------------------------------
  # 事後効果（alpha_hat 固定、rho の不確実性のみ）
  #   y_cf,t = alpha_hat' * (I - rho*(W + w alpha_hat'))^{-1}
  #            * [ (I - rho W) Yc_post,t - rho w Y0_post,t ]
  #--------------------------------------
  IN <- diag(N)
  w_l2 <- as.numeric(w)
  wn <- sqrt(sum(w_l2 * w_l2))
  if (wn > 0) {
    w_l2 <- w_l2 / wn
  }

  cf_one_rho <- function(rho) {
    Ainv <- solve(IN - rho * (W + w_l2 %*% t(alpha_hat)))
    B <- (IN - rho * W)
    ycf <- numeric(T1)
    for (t in seq_len(T1)) {
      tmp <- Ainv %*% (B %*% Yc_post[t, ] - rho * w_l2 * Y0_post[t])
      ycf[t] <- as.numeric(crossprod(alpha_hat, tmp))
    }
    ycf
  }

  spill_one_rho <- function(rho) {
    # Y_cf = (I - rho*W - rho*w*a')^{-1} * [ (I - rho*W)Yc - rho*w*Y0 ]
    
    # M_obs_inv = (I - rho*W - rho*w*a')^{-1}
    M_obs_inv <- solve(IN - rho * (W + w_l2 %*% t(alpha_hat)))
    
    # B = (I - rho*W)
    B <- (IN - rho * W)
    
    Yc_post_cf <- matrix(NA_real_, nrow = T1, ncol = N)
    for (t in seq_len(T1)) {
      # Yc_post_cf[t, ] = M_obs_inv %*% ( B %*% Yc_post[t, ] - rho * w_l2 * Y0_post[t] )
      
      # (cf_one_rho と同じ計算)
      tmp <- M_obs_inv %*% (B %*% Yc_post[t, ] - rho * w_l2 * Y0_post[t])
      Yc_post_cf[t, ] <- as.numeric(tmp)
    }
    
    # Spillover = Y_obs - Y_cf
    spill_effect <- Yc_post - Yc_post_cf
    spill_effect
  }

  # 点推定（rho_hat）
  ycf_point <- cf_one_rho(rho_hat)
  te_point <- as.numeric(Y0_post - ycf_point)
  ate_point <- mean(te_point)

  spill_draws_list <- lapply(rho_draws, spill_one_rho)
  # (T1 x N x M) の配列に変換
  spill_draws_array <- array(
    unlist(spill_draws_list), 
    dim = c(T1, N, length(rho_draws))
  )

  spill_mean_matrix <- apply(spill_draws_array, c(1, 2), mean, na.rm = TRUE)
  colnames(spill_mean_matrix) <- colnames(Yc_post) # ユニット名を付与

  # ATE の 95% CI（rho のみ回して近似）
  ate_draws <- vapply(
    rho_draws,
    function(r) {
      mean(Y0_post - cf_one_rho(r))
    },
    numeric(1)
  )
  ate_ci95 <- stats::quantile(
    ate_draws,
    c(0.025, 0.975),
    names = FALSE,
    type = 8
  )

  # 必要最小の効果要約を返す
  eff <- list(
    te_point = te_point, # T1-vector
    ate_point = ate_point, # scalar
    ate_ci95 = ate_ci95, # length-2
    spill = spill_mean_matrix
  )

  # 出力を従来の形に合わせて構築
  inputs <- list(
    Y0_pre = Y0_pre,
    Y0_post = Y0_post,
    Yc_pre = Yc_pre,
    Yc_post = Yc_post,
    times_pre = prep$times_pre,
    times_post = prep$times_post,
    units = prep$units,
    w = as.matrix(w),
    W = as.matrix(W)
  )

  structure(
    list(
      alpha_draws = alpha_draws, # Step1 の M x N
      rho_draws = rho_draws, # Step2 の M
      alpha_hat = alpha_hat, # colMeans(alpha_draws)
      rho_hat = rho_hat, # mean(rho_draws)
      effects = eff, # alpha_hat 固定・rho のみ反映
      inputs = inputs,
      sar = sar, # Step2 の詳細（rho/sigma2/Lambda/F など）
      T0 = T0
    ),
    class = "scspill"
  )
}

row_normalize <- function(W, tol = 1e-12, zero_policy = c("keep", "uniform0")) {
  zero_policy <- match.arg(zero_policy)
  stopifnot(is.matrix(W), is.numeric(W))
  W <- as.matrix(W)
  # 対角を0に（自己重みは使わない前提）
  diag(W) <- 0
  # 負の値があれば警告
  if (any(W < -tol, na.rm = TRUE)) {
    warning("W has negative entries.")
  }
  # NAは0として扱う
  W[is.na(W)] <- 0
  rs <- rowSums(W)
  # 0除算回避
  nz <- rs > tol
  W[nz, ] <- W[nz, , drop = FALSE] / rs[nz]
  # ゼロ行の扱い
  if (any(!nz)) {
    if (zero_policy == "uniform0") {
      # ゼロ行を一様(=0)のまま（既に0なので何もしない）
      # 何もしない
      NULL
    } else {
      # keep: 何もしない（同じ）
      NULL
    }
  }
  W
}

