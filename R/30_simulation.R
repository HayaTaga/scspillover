# 必要：quadprog（SCMの凸最適化）
if (!requireNamespace("quadprog", quietly = TRUE)) {
  install.packages("quadprog")
}

# ----------------------------
# 1回分のシミュレーション推定（SCM / BSCM / Proposed）
# ----------------------------
run_one_sim <- function(
  dgp,
  M_alpha = 2000, # α のMCMC反復
  burn_alpha = 1000, # α のバーンイン
  M_sar = 4000, # ρ のMCMC反復
  burn_sar = 2000, # ρ のバーンイン
  step_rho = 0.02, # ρ のRW幅（C++側と整合）
  add_pred_noise = TRUE # BSCMのポスト予測に pre 残差ノイズを足すか
) {
  Y0_pre <- dgp$data$Y0_pre
  Yc_pre <- dgp$data$Yc_pre
  Y0_post <- dgp$data$Y0_post
  Yc_post <- dgp$data$Yc_post
  W <- dgp$W
  w <- as.numeric(dgp$w)
  tau_true <- dgp$truth$tau_post
  T1 <- length(tau_true)
  N <- ncol(Yc_pre)

  ## =======================
  ## A) 古典的 SCM（preでα学習）
  ## =======================
  {
    y <- as.numeric(Y0_pre)
    X <- as.matrix(Yc_pre) # T0 x N
    ridge <- 1e-8
    D <- crossprod(X) + ridge * diag(N)
    d <- crossprod(X, y)
    Amat <- cbind(rep(1, N), diag(N)) # 制約: 1'a=1, a>=0
    bvec <- c(1, rep(0, N))
    qp <- quadprog::solve.QP(
      Dmat = D,
      dvec = d,
      Amat = Amat,
      bvec = bvec,
      meq = 1
    )
    a_scm <- drop(qp$solution)

    ycf_pre <- as.numeric(Yc_pre %*% a_scm)
    ycf_post <- as.numeric(Yc_post %*% a_scm)
    resid_pre <- Y0_pre - ycf_pre
    sd_pre <- stats::sd(resid_pre)

    eff_scm <- Y0_post - ycf_post
    lo_scm <- eff_scm - 1.96 * sd_pre
    hi_scm <- eff_scm + 1.96 * sd_pre
    covered_scm <- as.integer(lo_scm <= tau_true & tau_true <= hi_scm)

    ATE_scm_mean <- mean(eff_scm)
    ATE_scm_lo <- ATE_scm_mean - 1.96 * sd_pre / sqrt(T1)
    ATE_scm_hi <- ATE_scm_mean + 1.96 * sd_pre / sqrt(T1)
    ATE_scm_cov <- as.integer(
      ATE_scm_lo <= mean(tau_true) & mean(tau_true) <= ATE_scm_hi
    )

    scm <- list(
      alpha = a_scm,
      per_time = data.frame(
        t = seq_len(T1),
        eff_mean = eff_scm,
        eff_lo = lo_scm,
        eff_hi = hi_scm,
        covered = covered_scm
      ),
      ATE = data.frame(
        mean = ATE_scm_mean,
        lo = ATE_scm_lo,
        hi = ATE_scm_hi,
        covered = ATE_scm_cov
      )
    )
  }

  ## =======================
  ## B) Bayesian SCM（αのみ; preでMCMC, postで事後予測）
  ## =======================
  {
    alpha_draws <- hs_alpha_gibbs_cpp(
      Y0_pre = Y0_pre,
      control_outcome_pre = Yc_pre,
      iteration = M_alpha,
      burn = burn_alpha,
      verbose = FALSE
    ) # (M_alpha - burn) x N
    M_a <- nrow(alpha_draws)

    # pre残差sdを各ドローで推定（予測ノイズオプション用）
    sd_pre_by_m <- if (add_pred_noise) {
      apply(alpha_draws, 1, function(a) {
        r <- Y0_pre - as.numeric(Yc_pre %*% as.numeric(a))
        stats::sd(r)
      })
    } else {
      rep(NA_real_, M_a)
    }

    eff_draws <- matrix(NA_real_, nrow = T1, ncol = M_a)
    for (m in 1:M_a) {
      a <- as.numeric(alpha_draws[m, ])
      mu_post <- as.numeric(Yc_post %*% a)
      if (add_pred_noise) {
        eff_draws[, m] <- Y0_post - (mu_post + rnorm(T1, 0, sd_pre_by_m[m]))
      } else {
        eff_draws[, m] <- Y0_post - mu_post
      }
    }
    eff_mean <- rowMeans(eff_draws)
    q_b <- t(apply(eff_draws, 1, stats::quantile, probs = c(0.025, 0.975)))
    covered_b <- as.integer(q_b[, 1] <= tau_true & tau_true <= q_b[, 2])

    ATE_b_means <- colMeans(eff_draws)
    ATE_b_lohi <- stats::quantile(ATE_b_means, c(0.025, 0.975))
    ATE_b_cov <- as.integer(
      ATE_b_lohi[1] <= mean(tau_true) & mean(tau_true) <= ATE_b_lohi[2]
    )

    bscm <- list(
      alpha_draws = alpha_draws,
      per_time = data.frame(
        t = seq_len(T1),
        eff_mean = eff_mean,
        eff_lo = q_b[, 1],
        eff_hi = q_b[, 2],
        covered = covered_b
      ),
      ATE = data.frame(
        mean = mean(ATE_b_means),
        lo = ATE_b_lohi[1],
        hi = ATE_b_lohi[2],
        covered = ATE_b_cov
      )
    )
  }

  ## =======================
  ## C) Proposed: SC-SPILL（αとρを併用）
  ##    αはBSCMのドロー、ρはSAR MCMCのドローを使用し、同時式で反事実を再構成
  ## =======================
  {
    # SAR（pre; K=0, p=0で十分）
    sar <- sar_full_sampler_cpp(
      Y0_pre = Y0_pre,
      Yc_pre = Yc_pre,
      Xc_pre_ = NULL, # K=0
      T0 = nrow(Yc_pre),
      N = ncol(Yc_pre),
      K = 0,
      p = 0,
      w = w,
      W = W,
      iteration = M_sar,
      burn = burn_sar,
      step_rho = step_rho,
      a0 = 1.0,
      b0 = 1.0,
      verbose = FALSE
    )
    rho_draws <- as.numeric(sar$rho) # (M_sar - burn_sar)
    M_s <- length(rho_draws)

    # α と ρ のドロー数を合わせる（短い方に合わせて先頭を使用）
    M <- min(nrow(bscm$alpha_draws), M_s)
    A <- as.matrix(bscm$alpha_draws[1:M, , drop = FALSE]) # M x N
    R <- rho_draws[1:M]

    IN <- diag(N)
    eff_p_draws <- matrix(NA_real_, nrow = T1, ncol = M)

    # 1時点ずつ： y_cf = a' Ainv { (I - r W) y_c - r w y0 }
    for (t in 1:T1) {
      yc <- as.numeric(Yc_post[t, ])
      y0 <- Y0_post[t]
      for (m in 1:M) {
        a <- as.numeric(A[m, ])
        r <- R[m]
        # 数値安定（行列が特異ならNA）
        Mmat <- tryCatch(IN - r * (w %*% t(a) + W), error = function(e) NULL)
        if (is.null(Mmat)) {
          eff_p_draws[t, m] <- NA_real_
          next
        }
        Binv <- tryCatch(solve(Mmat), error = function(e) NULL)
        if (is.null(Binv)) {
          eff_p_draws[t, m] <- NA_real_
          next
        }
        tmp <- Binv %*% ((IN - r * W) %*% yc - r * w * y0)
        y_cf <- sum(a * as.numeric(tmp))
        eff_p_draws[t, m] <- y0 - y_cf
      }
    }

    eff_p_mean <- rowMeans(eff_p_draws, na.rm = TRUE)
    q_p <- t(apply(
      eff_p_draws,
      1,
      stats::quantile,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    ))
    covered_p <- as.integer(q_p[, 1] <= tau_true & tau_true <= q_p[, 2])

    ATE_p_means <- colMeans(eff_p_draws, na.rm = TRUE)
    ATE_p_lohi <- stats::quantile(ATE_p_means, c(0.025, 0.975), na.rm = TRUE)
    ATE_p_cov <- as.integer(
      ATE_p_lohi[1] <= mean(tau_true) & mean(tau_true) <= ATE_p_lohi[2]
    )

    prop <- list(
      alpha_rho_M = M,
      per_time = data.frame(
        t = seq_len(T1),
        eff_mean = eff_p_mean,
        eff_lo = q_p[, 1],
        eff_hi = q_p[, 2],
        covered = covered_p
      ),
      ATE = data.frame(
        mean = mean(ATE_p_means, na.rm = TRUE),
        lo = ATE_p_lohi[1],
        hi = ATE_p_lohi[2],
        covered = ATE_p_cov
      )
    )
  }

  list(
    truth = list(alpha = dgp$truth$alpha, tau = tau_true, ATE = mean(tau_true)),
    scm = scm,
    bscm = bscm,
    prop = prop
  )
}

# ----------------------------
# 多回シミュレーションの集計（3方式）
# ----------------------------
summarize_many <- function(results) {
  stopifnot(length(results) >= 1)
  tau_true <- results[[1]]$truth$tau
  T1 <- length(tau_true)
  ATE_true <- results[[1]]$truth$ATE

  # 期別の推定値と被覆
  mat_extract <- function(path) {
    do.call(cbind, lapply(results, function(r) r[[path]]$per_time$eff_mean))
  }
  cov_extract <- function(path) {
    do.call(cbind, lapply(results, function(r) r[[path]]$per_time$covered))
  }

  scm_eff <- mat_extract("scm")
  bscm_eff <- mat_extract("bscm")
  prop_eff <- mat_extract("prop")

  scm_cov <- cov_extract("scm")
  bscm_cov <- cov_extract("bscm")
  prop_cov <- cov_extract("prop")

  # 期別の平均バイアス/RMSE/カバレッジ
  bias_rmse_cov <- function(eff_mat, cov_mat) {
    bias_t <- rowMeans(eff_mat - tau_true)
    rmse_t <- sqrt(rowMeans((eff_mat - tau_true)^2))
    cover_t <- rowMeans(cov_mat)
    list(bias_t = bias_t, rmse_t = rmse_t, cover_t = cover_t)
  }
  s_scm <- bias_rmse_cov(scm_eff, scm_cov)
  s_bscm <- bias_rmse_cov(bscm_eff, bscm_cov)
  s_prop <- bias_rmse_cov(prop_eff, prop_cov)

  # ATE 集計
  ATE_vec <- function(path) sapply(results, function(r) r[[path]]$ATE$mean)
  ATE_cov <- function(path) sapply(results, function(r) r[[path]]$ATE$covered)

  scm_ATE <- ATE_vec("scm")
  scm_ATE_cov <- ATE_cov("scm")
  bscm_ATE <- ATE_vec("bscm")
  bscm_ATE_cov <- ATE_cov("bscm")
  prop_ATE <- ATE_vec("prop")
  prop_ATE_cov <- ATE_cov("prop")

  summary_methods <- rbind(
    data.frame(
      method = "SCM",
      ATE_truth = ATE_true,
      ATE_mean_est = mean(scm_ATE),
      ATE_bias = mean(scm_ATE - ATE_true),
      ATE_rmse = sqrt(mean((scm_ATE - ATE_true)^2)),
      ATE_cover95 = mean(scm_ATE_cov),
      per_time_avg_bias = mean(s_scm$bias_t),
      per_time_avg_rmse = mean(s_scm$rmse_t),
      per_time_avg_cover95 = mean(s_scm$cover_t),
      row.names = NULL
    ),
    data.frame(
      method = "Bayesian SCM",
      ATE_truth = ATE_true,
      ATE_mean_est = mean(bscm_ATE),
      ATE_bias = mean(bscm_ATE - ATE_true),
      ATE_rmse = sqrt(mean((bscm_ATE - ATE_true)^2)),
      ATE_cover95 = mean(bscm_ATE_cov),
      per_time_avg_bias = mean(s_bscm$bias_t),
      per_time_avg_rmse = mean(s_bscm$rmse_t),
      per_time_avg_cover95 = mean(s_bscm$cover_t),
      row.names = NULL
    ),
    data.frame(
      method = "Proposed (SC-SPILL)",
      ATE_truth = ATE_true,
      ATE_mean_est = mean(prop_ATE),
      ATE_bias = mean(prop_ATE - ATE_true),
      ATE_rmse = sqrt(mean((prop_ATE - ATE_true)^2)),
      ATE_cover95 = mean(prop_ATE_cov),
      per_time_avg_bias = mean(s_prop$bias_t),
      per_time_avg_rmse = mean(s_prop$rmse_t),
      per_time_avg_cover95 = mean(s_prop$cover_t),
      row.names = NULL
    )
  )

  list(
    summary_methods = summary_methods,
    per_time = list(
      scm = data.frame(
        t = seq_len(T1),
        bias = s_scm$bias_t,
        rmse = s_scm$rmse_t,
        cover95 = s_scm$cover_t
      ),
      bscm = data.frame(
        t = seq_len(T1),
        bias = s_bscm$bias_t,
        rmse = s_bscm$rmse_t,
        cover95 = s_bscm$cover_t
      ),
      prop = data.frame(
        t = seq_len(T1),
        bias = s_prop$bias_t,
        rmse = s_prop$rmse_t,
        cover95 = s_prop$cover_t
      )
    )
  )
}
