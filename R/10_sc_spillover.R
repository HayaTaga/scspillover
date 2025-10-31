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
  step_rho = 0.05
) {
  stopifnot(is.data.frame(data))
  if (!is.character(y)) {
    y <- as.character(substitute(y))
  }

  # 必須列チェック
  required_cols <- c(unit_col, time_col, y, treatment_dummy)
  if (!all(required_cols %in% names(data))) {
    stop(sprintf(
      "`data` に列 %s が必要です。",
      paste0(required_cols, collapse = ", ")
    ))
  }
  set.seed(seed)

  # 介入開始期から T0 を決める
  times <- sort(unique(data[[time_col]]))
  treat_series <- data[
    data[[unit_col]] == treated_unit,
    c(time_col, treatment_dummy)
  ]
  t_start <- min(treat_series[[time_col]][treat_series[[treatment_dummy]] == 1])
  T0 <- if (is.null(T0)) sum(times < t_start) else as.integer(T0)

  # 整形 (Y0_pre, Yc_pre, ...) を得る
  prep <- scspill_prep_X(
    data,
    treated_unit = treated_unit,
    T0 = T0,
    X = X,
    y_col = y,
    unit_col = unit_col,
    time_col = time_col
  )

  joint_norm <- normalize_joint_wW(W, w)
  W_use <- joint_norm$W_new
  w_use <- joint_norm$w_new

  Y0_pre <- prep$Y0_pre
  Yc_pre <- prep$Yc_pre
  Y0_post <- prep$Y0_post
  Yc_post <- prep$Yc_post
  Xc_pre <- prep$Xc_pre

  N <- ncol(Yc_pre)
  T1 <- nrow(Yc_post)

  #--------------------------------------
  # Step 1: BSCM による alpha の推定
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
  # Step 2: α̂ 固定で rho 等を推定
  #--------------------------------------
  if (verbose) {
    message("[Step 2] Sampling rho (and others) with fixed alpha_hat...")
  }

  # X の整形
  K <- 0L
  Xvec <- NULL
  if (!is.null(Xc_pre)) {
    if (is.array(Xc_pre) && length(dim(Xc_pre)) == 3L) {
      K <- dim(Xc_pre)[3]
      stopifnot(dim(Xc_pre)[1] == T0, dim(Xc_pre)[2] == N)
      Xvec <- as.numeric(aperm(Xc_pre, c(1, 2, 3)))
    } else {
      Xvec <- as.numeric(Xc_pre)
      Ktmp <- length(Xvec) / (T0 * N)
      if (abs(Ktmp - round(Ktmp)) > 1e-8) {
        stop("Xc_pre の次元が (T0*N*K) に整合しません。")
      }
      K <- as.integer(round(Ktmp))
    }
  }

  w_l2 <- as.numeric(w_use)

  sar <- sar_full_sampler_cpp_step2(
    Yc_pre = Yc_pre,
    alpha_hat_in = alpha_hat,
    Xc_pre_ = if (!is.null(Xvec)) Xvec else R_NilValue,
    T0 = T0,
    N = N,
    K = K,
    p = as.integer(p_factors),
    w_in = w_l2,
    W = as.matrix(W_use),
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
  #--------------------------------------
  IN <- diag(N)

  cf_one_rho <- function(rho) {
    Ainv <- solve(IN - rho * (W_use + w_l2 %*% t(alpha_hat)))
    B <- (IN - rho * W_use)
    ycf <- numeric(T1)
    for (t in seq_len(T1)) {
      tmp <- Ainv %*% (B %*% Yc_post[t, ] - rho * w_l2 * Y0_post[t])
      ycf[t] <- as.numeric(crossprod(alpha_hat, tmp))
    }
    ycf
  }

  spill_one_rho <- function(rho) {
    M_obs_inv <- solve(IN - rho * (W_use + w_l2 %*% t(alpha_hat)))
    B <- (IN - rho * W_use)

    Yc_post_cf <- matrix(NA_real_, nrow = T1, ncol = N)
    for (t in seq_len(T1)) {
      tmp <- M_obs_inv %*% (B %*% Yc_post[t, ] - rho * w_l2 * Y0_post[t])
      Yc_post_cf[t, ] <- as.numeric(tmp)
    }

    spill_effect <- Yc_post - Yc_post_cf
    spill_effect
  }

  ycf_point <- cf_one_rho(rho_hat)
  te_point <- as.numeric(Y0_post - ycf_point)
  ate_point <- mean(te_point)

  spill_draws_list <- lapply(rho_draws, spill_one_rho)
  spill_draws_array <- array(
    unlist(spill_draws_list),
    dim = c(T1, N, length(rho_draws))
  )
  spill_mean_matrix <- apply(spill_draws_array, c(1, 2), mean, na.rm = TRUE)
  colnames(spill_mean_matrix) <- colnames(Yc_post)

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
    names = FALSE
  )

  eff <- list(
    te_point = te_point,
    ate_point = ate_point,
    ate_ci95 = ate_ci95,
    spill = spill_mean_matrix
  )

  inputs <- list(
    Y0_pre = Y0_pre,
    Y0_post = Y0_post,
    Yc_pre = Yc_pre,
    Yc_post = Yc_post,
    times_pre = prep$times_pre,
    times_post = prep$times_post,
    units = prep$units,
    w = as.matrix(w_use),
    W = as.matrix(W_use)
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

#' @export
normalize_joint_wW <- function(W, w) {
  stopifnot(is.matrix(W))
  stopifnot(length(w) == nrow(W), ncol(W) == nrow(W))

  N <- nrow(W)
  joint <- cbind(w, W)

  rs <- rowSums(joint)
  rs[rs == 0] <- 1

  joint_norm <- joint / rs

  w_new <- joint_norm[, 1, drop = TRUE]
  W_new <- joint_norm[, -1, drop = FALSE]

  list(W_new = W_new, w_new = w_new)
}
