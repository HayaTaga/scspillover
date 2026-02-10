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

  required_cols <- c(unit_col, time_col, y, treatment_dummy)
  if (!all(required_cols %in% names(data))) {
    stop(sprintf(
      "Required columns missing in data: %s",
      paste0(required_cols, collapse = ", ")
    ))
  }
  set.seed(seed)

  times <- sort(unique(data[[time_col]]))
  treat_series <- data[
    data[[unit_col]] == treated_unit,
    c(time_col, treatment_dummy)
  ]
  t_start <- min(treat_series[[time_col]][treat_series[[treatment_dummy]] == 1])
  T0 <- if (is.null(T0)) sum(times < t_start) else as.integer(T0)

  prep <- scspill_prep_X(
    data,
    treated_unit = treated_unit,
    T0 = T0,
    X = X,
    y_col = y,
    unit_col = unit_col,
    time_col = time_col
  )

  W_use <- row_normalize(W)
  w_use <- as.numeric(w)
  wsum <- sum(w_use)
  if (is.finite(wsum) && wsum > 1e-12) {
    w_use <- w_use / wsum
  }
  Y0_pre <- prep$Y0_pre
  Yc_pre <- prep$Yc_pre
  Y0_post <- prep$Y0_post
  Yc_post <- prep$Yc_post
  Xc_pre <- prep$Xc_pre

  N <- ncol(Yc_pre)
  T1 <- nrow(Yc_post)

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

  if (verbose) {
    message("[Step 2] Sampling rho (and others) with fixed alpha_hat...")
  }

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
        stop("Xc_pre dimensions do not match (T0*N*K).")
      }
      K <- as.integer(round(Ktmp))
    }
  }

  sar <- sar_full_sampler_cpp_step2(
    Yc_pre = Yc_pre,
    alpha_hat_in = alpha_hat,
    Xc_pre_ = if (!is.null(Xvec)) Xvec else R_NilValue,
    T0 = T0,
    N = N,
    K = K,
    p = as.integer(p_factors),
    w = w_use,
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

  IN <- diag(N)

  cf_one_rho <- function(rho) {
    Ainv <- solve(IN - rho * (W_use + w_use %*% t(alpha_hat)))
    B <- (IN - rho * W_use)
    ycf <- numeric(T1)
    for (t in seq_len(T1)) {
      tmp <- Ainv %*% (B %*% Yc_post[t, ] - rho * w_use * Y0_post[t])
      ycf[t] <- as.numeric(crossprod(alpha_hat, tmp))
    }
    ycf
  }

  spill_one_rho <- function(rho) {
    M_obs_inv <- solve(IN - rho * (W_use + w_use %*% t(alpha_hat)))
    B <- (IN - rho * W_use)

    Yc_pre_cf <- matrix(NA_real_, nrow = T0, ncol = N)
    for (t in seq_len(T0)) {
      tmp <- M_obs_inv %*% (B %*% Yc_pre[t, ] - rho * w_use * Y0_pre[t])
      Yc_pre_cf[t, ] <- as.numeric(tmp)
    }
    spill_pre <- Yc_pre - Yc_pre_cf

    Yc_post_cf <- matrix(NA_real_, nrow = T1, ncol = N)
    for (t in seq_len(T1)) {
      tmp <- M_obs_inv %*% (B %*% Yc_post[t, ] - rho * w_use * Y0_post[t])
      Yc_post_cf[t, ] <- as.numeric(tmp)
    }
    spill_post <- Yc_post - Yc_post_cf

    rbind(spill_pre, spill_post)
  }

  ycf_point <- cf_one_rho(rho_hat)
  te_point <- as.numeric(Y0_post - ycf_point)
  ate_point <- mean(te_point)

  spill_draws_list <- lapply(rho_draws, spill_one_rho)
  spill_draws_array <- array(
    unlist(spill_draws_list),
    dim = c(T0 + T1, N, length(rho_draws))
  )

  spill_mean_matrix <- apply(spill_draws_array, c(1, 2), mean, na.rm = TRUE)
  colnames(spill_mean_matrix) <- colnames(Yc_post)
  rownames(spill_mean_matrix) <- c(prep$times_pre, prep$times_post)

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

#' Row-normalize a spatial weights matrix W
#'
#' Ensures that the sum of each row is 1.
#' Sets the diagonal to 0 and handles rows that sum to 0.
#'
#' @param W A numeric matrix.
#' @param tol Tolerance for checking if a row sum is zero.
#' @param zero_policy How to handle rows that sum to zero (or are close to it).
#'   "keep" (default): leaves the row as all zeros.
#'   "uniform": (not implemented here, but common) sets to 1/N.
#' @return A row-normalized matrix.
#'
row_normalize <- function(W, tol = 1e-12, zero_policy = c("keep")) {
  zero_policy <- match.arg(zero_policy)

  if (!is.matrix(W) || !is.numeric(W)) {
    stop("W must be a numeric matrix.")
  }

  # Ensure diagonal is zero (no self-loops)
  diag(W) <- 0

  # Calculate row sums
  rs <- rowSums(W, na.rm = TRUE)

  # Find rows that are not zero (or very close to it)
  nz <- rs > tol

  # Normalize non-zero rows
  if (any(nz)) {
    W[nz, ] <- W[nz, , drop = FALSE] / rs[nz]
  }

  # Handle zero-sum rows (if any)
  if (any(!nz)) {
    if (zero_policy == "keep") {
      # Do nothing, leave the row as all zeros
      W[!nz, ] <- 0 # Ensure it's clean
    }
  }

  W
}
