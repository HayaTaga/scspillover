#' Synthetic Control with Spillovers (Bayesian, compact)
#'
#' @param data long panel with columns unit, time, y (+ optional covariates)
#' @param treated_unit treated unit name
#' @param T0 number of pre-treatment periods
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
     treatment_dummy
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
     T0 <- sum(times < t_start)

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

     w <- row_normalize(w)
     W <- row_normalize(W)

     Y0_pre <- prep$Y0_pre
     Yc_pre <- prep$Yc_pre
     Y0_post <- prep$Y0_post
     Yc_post <- prep$Yc_post
     Xc_pre <- prep$Xc_pre
     N <- ncol(Yc_pre)
     T1 <- nrow(Yc_post)

     # (i) α推定（horseshoe）
     # hs <- hs_alpha_gibbs(
     #      y = Y0_pre,
     #      X = Yc_pre,
     #      M = M,
     #      burn = burn,
     #      verbose = verbose
     # )
     hs <- hs_alpha_gibbs(
          y = Y0_pre,
          X = Yc_pre,
          M = M,
          burn = burn,
          verbose = verbose
     )
     alpha_draws <- hs$alpha
     alpha_hat <- colMeans(alpha_draws)

     # (ii) SAR推定（ρ, β）
     # sar <- sar_gibbs_sampler(
     #      Y0_pre = Y0_pre,
     #      Yc_pre = Yc_pre,
     #      Xc_pre = Xc_pre,
     #      W = W,
     #      w = w,
     #      M = M,
     #      burn = burn,
     #      verbose = verbose,
     #      p_factors = p_factors
     # )
     sar <- sar_gibbs_sampler(
          Y0_pre = Y0_pre,
          Yc_pre = Yc_pre,
          Xc_pre = Xc_pre,
          W = W,
          w = w,
          M = M,
          burn = burn,
          verbose = verbose,
          p_factors = p_factors
     )
     rho_draws <- sar$rho
     rho_hat <- mean(rho_draws)

     # (iii) 効果算出（式(5)(6)）
     eff <- posterior_effects(
          Y0_post = Y0_post,
          Yc_post = Yc_post,
          alpha_draws = alpha_draws,
          rho_draws = rho_draws,
          w = w,
          W = W
     )

     # structure(
     #      list(
     #           alpha_draws = alpha_draws,
     #           alpha_hat = alpha_hat,
     #           rho_draws = rho_draws,
     #           rho_hat = rho_hat,
     #           effects = eff,
     #           sar = sar,
     #           inputs = list(
     #                units = prep$units,
     #                times = prep$times,
     #                T0 = T0,
     #                N = N,
     #                T1 = T1,
     #                X_used = !is.null(Xc_pre),
     #                p_factors = p_factors,
     #                y_col = y,
     #                unit_col = unit_col,
     #                time_col = time_col
     #           )
     #      ),
     #      class = "scspill"
     # )
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

     structure(list(
          alpha_draws = alpha_draws,
          rho_draws = rho_draws,
          alpha_hat = alpha_hat,
          rho_hat = rho_hat,
          effects = eff,
          inputs = inputs,
          sar = sar,
          T0 = T0,
          hs=hs,
          sar=sar
          ),
          class = "scspill"
     )
}


row_normalize <- function(x, margin = 1, tol = 1e-12) {
     if (is.null(x)) {
          return(NULL)
     }
     if (is.vector(x)) {
          s <- sum(x, na.rm = TRUE)
          if (abs(s) < tol) {
               return(rep(0, length(x)))
          }
          return(x / s)
     }
     if (is.matrix(x)) {
          if (margin == 1) {
               row_sums <- rowSums(x, na.rm = TRUE)
               row_sums[row_sums < tol] <- 1
               return(x / row_sums)
          } else if (margin == 2) {
               col_sums <- colSums(x, na.rm = TRUE)
               col_sums[col_sums < tol] <- 1
               return(t(t(x) / col_sums))
          } else {
               stop("margin must be 1 (rows) or 2 (columns)")
          }
     }
     stop("x must be numeric vector or matrix")
}
