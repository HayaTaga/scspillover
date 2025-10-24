#' Synthetic Control with Spillovers (Bayesian, joint MCMC)
#'
#' @param data long panel with columns unit, time, y (+ optional covariates)
#' @param treated_unit treated unit name
#' @param T0 number of pre-treatment periods (if NULL, inferred from treatment start)
#' @param w length-N vector (controls->treated)
#' @param W N x N row-normalized spatial weight among controls
#' @param X NULL, character vector of covariate colnames in `data`, or matrix with nrow = nrow(data)
#' @param p_factors number of latent factors (>=0; 0 で因子層無効)
#' @param M posterior draws; @param burn burn-in; @param seed RNG; @param verbose logical
#' @export
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
     step_alpha = 0.05
) {
     stopifnot(is.data.frame(data))
     if (!is.character(y)) {
          y <- as.character(substitute(y))
     }

     # 存在チェック
     required_cols <- c(unit_col, time_col, y, treatment_dummy)
     if (!all(required_cols %in% names(data))) {
          stop(sprintf(
               "`data` に列 %s が必要です。",
               paste0(required_cols, collapse = ", ")
          ))
     }

     set.seed(seed)

     times <- sort(unique(data[[time_col]]))
     treat_series <- data[
          data[[unit_col]] == treated_unit,
          c(time_col, treatment_dummy)
     ]
     t_start <- min(treat_series[[time_col]][
          treat_series[[treatment_dummy]] == 1
     ])
     T0 <- if (is.null(T0)) sum(times < t_start) else as.integer(T0)

     # 前処理（列名を指定）
     prep <- scspill_prep_X(
          data,
          treated_unit = treated_unit,
          T0 = T0,
          X = X,
          y_col = y,
          unit_col = unit_col,
          time_col = time_col
     )

     # 正規化（ゼロ割り回避ロバスト）
     w <- row_normalize(w)
     W <- row_normalize(W)

     Y0_pre <- prep$Y0_pre
     Yc_pre <- prep$Yc_pre
     Y0_post <- prep$Y0_post
     Yc_post <- prep$Yc_post
     Xc_pre <- prep$Xc_pre
     N <- ncol(Yc_pre)
     T1 <- nrow(Yc_post)

     # ================================
     #  共同サンプラー（rho, alpha を同時推定）
     # ================================
     sar <- sar_gibbs_sampler(
          Y0_pre = Y0_pre,
          Yc_pre = Yc_pre,
          Xc_pre = Xc_pre,
          W = W,
          w = w,
          M = M,
          burn = burn,
          verbose = verbose,
          p_factors = p_factors,
          step_rho = step_rho,
          step_alpha = step_alpha
     )

     rho_draws <- sar$rho # length M
     alpha_draws <- sar$alpha # M x N
     rho_hat <- mean(rho_draws)
     alpha_hat <- colMeans(alpha_draws)

     # ================================
     #  事後効果（式 (5)(6) に基づく）
     # ================================
     eff <- posterior_effects(
          Y0_post = Y0_post,
          Yc_post = Yc_post,
          alpha_draws = alpha_draws,
          rho_draws = rho_draws,
          w = w,
          W = W
     )

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
               alpha_draws = alpha_draws,
               rho_draws = rho_draws,
               alpha_hat = alpha_hat,
               rho_hat = rho_hat,
               effects = eff,
               inputs = inputs,
               sar = sar,
               T0 = T0
          ),
          class = "scspill"
     )
}

row_normalize <- function(W, tol = 1e-12, zero_policy = c("keep","uniform0")) {
  zero_policy <- match.arg(zero_policy)
  stopifnot(is.matrix(W), is.numeric(W))
  W <- as.matrix(W)
  # 対角を0に（自己重みは使わない前提）
  diag(W) <- 0
  # 負の値があれば警告
  if (any(W < -tol, na.rm = TRUE)) warning("W has negative entries.")
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