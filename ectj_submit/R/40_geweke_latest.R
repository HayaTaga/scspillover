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


#' @keywords internal
compute_bnd <- function(W_use, w_use, alpha_hat_scaled, c_stability = 0.95) {
  A <- W_use + w_use %*% t(alpha_hat_scaled)
  ev <- eigen(A, symmetric = FALSE, only.values = TRUE)$values
  maxabs <- max(Mod(ev))
  if (!is.finite(maxabs) || maxabs < 1e-12) {
    maxabs <- 1e-12
  }
  c_stability / maxabs
}

# --- g(θ, y) ---
#' @keywords internal
default_g_fn <- function(theta, Yc, Y0_pre, W_use, w_use) {
  yc_vec <- as.numeric(Yc)
  N <- ncol(Yc)
  T0 <- nrow(Yc)
  wyc <- as.numeric(Yc %*% as.numeric(w_use))
  spatial_q <- sum(diag(Yc %*% W_use %*% t(Yc))) / (N * T0)
  corr_y0_wyc <- if (stats::sd(Y0_pre) > 0 && stats::sd(wyc) > 0) {
    cor(Y0_pre, wyc)
  } else {
    NA_real_
  }

  beta_mean <- if (length(theta$beta) > 0) mean(theta$beta) else NA_real_
  Eta_mean <- if (length(theta$Eta) > 0) mean(theta$Eta) else NA_real_
  Gamma_mean <- if (length(theta$Gamma) > 0) mean(theta$Gamma) else NA_real_

  c(
    rho = unname(theta$rho),
    log_sigma2 = log(pmax(theta$sigma2, 1e-12)),
    yc_mean = mean(yc_vec),
    log_yc_var = log(pmax(stats::var(yc_vec), 1e-12)),
    spatial_quadratic = spatial_q,
    corr_y0_wyc = corr_y0_wyc,
    beta_mean = beta_mean,
    Eta_mean = Eta_mean,
    Gamma_mean = Gamma_mean
  )
}

#' @keywords internal
draw_initial_state <- function(
    T0,
    N,
    K,
    p,
    a0,
    b0,
    W_use,
    w_use,
    alpha_hat_scaled) {
  bnd <- compute_bnd(W_use, w_use, alpha_hat_scaled)
  rho0 <- stats::runif(1, -bnd, bnd)
  sigma2_0 <- 1 / stats::rgamma(1, shape = a0, rate = b0)
  beta0 <- if (K > 0) stats::rnorm(K, 0, 1) else numeric()
  Eta0 <- if (p > 0) {
    matrix(stats::rnorm(N * p, 0, 1), N, p)
  } else {
    matrix(0, N, 0)
  }
  Gamma0 <- if (p > 0) {
    matrix(stats::rnorm(p * T0, 0, 1), p, T0)
  } else {
    matrix(0, 0, T0)
  }
  list(
    rho = as.numeric(rho0),
    sigma2 = as.numeric(sigma2_0),
    beta = beta0,
    Eta = Eta0,
    Gamma = Gamma0
  )
}

#' @keywords internal
var_mcmc_batchmeans <- function(x, b = NULL) {
  x <- x[is.finite(x)]
  M <- length(x)
  if (is.null(b)) {
    b <- max(2L, floor(sqrt(M)))
  }
  a <- floor(M / b)
  if (a < 2L) {
    vx <- stats::var(x)
    return(list(var_mean = vx / M))
  }
  x2 <- x[seq_len(a * b)]
  bm <- colMeans(matrix(x2, nrow = b, ncol = a))
  tau2 <- b * stats::var(bm)
  list(var_mean = as.numeric(tau2) / (a * b))
}

#' @keywords internal
geweke_jdt_full <- function(
    Y0_pre,
    Yc_pre_like_dims, # c(T0, N)
    W,
    w,
    alpha_hat_scaled, # N
    Xc_pre = NULL, # T0 x N x K array or NULL
    p = 0L,
    M1 = 20000L,
    M2 = 20000L,
    burn_in = 5000L,
    a0 = 1.0,
    b0 = 1.0,
    step_rho = 0.05,
    g_fn = default_g_fn,
    batch_size = NULL,
    verbose = TRUE,
    rho_support = NULL) {
  stopifnot(is.numeric(Y0_pre))
  T0 <- length(Y0_pre)
  N <- Yc_pre_like_dims[2]
  K <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  W_use <- row_normalize(W)
  w_use <- as.numeric(w)
  wsum <- sum(w_use)
  if (is.finite(wsum) && wsum > 1e-12) {
    w_use <- w_use / wsum
  }

  spectral_bound <- function(W) {
    ev <- eigen(W, symmetric = FALSE, only.values = TRUE)$values
    r <- max(Mod(ev))
    if (!is.finite(r) || r <= 0) {
      return(0.95)
    }
    0.95 / r
  }
  if (is.null(rho_support)) {
    bnd <- spectral_bound(W_use)
    rho_lo <- -bnd
    rho_hi <- bnd
  } else {
    stopifnot(length(rho_support) == 2L, rho_support[1] < rho_support[2])
    rho_lo <- as.numeric(rho_support[1])
    rho_hi <- as.numeric(rho_support[2])
  }

  # ---------- MC (iid) side ----------
  if (verbose) {
    message("[JDT] MC side (iid) ...")
  }
  g_iid_mat <- NULL
  for (m in seq_len(M1)) {
    st <- draw_initial_state(
      T0,
      N,
      K,
      p,
      a0,
      b0,
      W_use,
      w_use,
      alpha_hat_scaled
    )
    Xc_used <- if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      st$rho,
      st$sigma2,
      Xc_used,
      if (K > 0) st$beta else numeric(),
      if (p > 0) st$Eta else matrix(0, N, 0),
      if (p > 0) st$Gamma else matrix(0, 0, T0)
    )
    g_val <- g_fn(st, Yc_draw, Y0_pre, W_use, w_use)
    if (is.null(g_iid_mat)) {
      g_iid_mat <- matrix(NA_real_, nrow = M1, ncol = length(g_val))
      colnames(g_iid_mat) <- names(g_val)
    }
    g_iid_mat[m, ] <- g_val
  }

  # ---------- SC (successive-conditional) side ----------
  if (verbose) {
    message("[JDT] SC side (successive-conditional) ...")
  }
  state <- draw_initial_state(
    T0,
    N,
    K,
    p,
    a0,
    b0,
    W_use,
    w_use,
    alpha_hat_scaled
  )

  # burn-in
  for (m in seq_len(burn_in)) {
    Xc_used <- if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      state$rho,
      state$sigma2,
      Xc_used,
      if (K > 0) state$beta else numeric(),
      if (p > 0) state$Eta else matrix(0, N, 0),
      if (p > 0) state$Gamma else matrix(0, 0, T0)
    )
    state <- scspill_one_step_cpp(
      Yc_data = Yc_draw,
      W_use = W_use,
      w_use = w_use,
      alpha_hat_scaled = alpha_hat_scaled,
      T0 = T0,
      N = N,
      Xc_pre = if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre,
      K = K,
      p = p,
      state_in = state,
      a0 = a0,
      b0 = b0,
      step_rho = step_rho,
      rho_lo = rho_lo,
      rho_hi = rho_hi
    )
  }

  # keep
  g_mcmc_mat <- matrix(NA_real_, nrow = M2, ncol = ncol(g_iid_mat))
  colnames(g_mcmc_mat) <- colnames(g_iid_mat)

  for (m in seq_len(M2)) {
    Xc_used <- if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre
    Yc_draw <- simulate_Yc_forward_cpp(
      T0,
      W_use,
      w_use,
      alpha_hat_scaled,
      state$rho,
      state$sigma2,
      Xc_used,
      if (K > 0) state$beta else numeric(),
      if (p > 0) state$Eta else matrix(0, N, 0),
      if (p > 0) state$Gamma else matrix(0, 0, T0)
    )
    state <- scspill_one_step_cpp(
      Yc_data = Yc_draw,
      W_use = W_use,
      w_use = w_use,
      alpha_hat_scaled = alpha_hat_scaled,
      T0 = T0,
      N = N,
      Xc_pre = if (is.null(Xc_pre)) array(0, c(T0, N, 0L)) else Xc_pre,
      K = K,
      p = p,
      state_in = state,
      a0 = a0,
      b0 = b0,
      step_rho = step_rho,
      rho_lo = rho_lo,
      rho_hi = rho_hi
    )
    g_mcmc_mat[m, ] <- g_fn(state, Yc_draw, Y0_pre, W_use, w_use)
  }

  mean_iid <- colMeans(g_iid_mat, na.rm = TRUE)
  mean_mcmc <- colMeans(g_mcmc_mat, na.rm = TRUE)

  n_iid <- colSums(is.finite(g_iid_mat))
  se_iid <- sqrt(
    apply(g_iid_mat, 2, stats::var, na.rm = TRUE) / pmax(n_iid, 1L)
  )

  if (is.null(batch_size)) {
    batch_size <- max(2L, floor(sqrt(M2)))
  }
  se_mcmc <- vapply(
    seq_len(ncol(g_mcmc_mat)),
    function(j) sqrt(var_mcmc_batchmeans(g_mcmc_mat[, j], batch_size)$var_mean),
    numeric(1L)
  )

  Z <- (mean_iid - mean_mcmc) / sqrt(se_iid^2 + se_mcmc^2)
  pval <- 2 * stats::pnorm(-abs(Z))

  summary <- data.frame(
    g = colnames(g_iid_mat),
    mean_iid = as.numeric(mean_iid),
    mean_mcmc = as.numeric(mean_mcmc),
    se_iid = as.numeric(se_iid),
    se_mcmc = as.numeric(se_mcmc),
    Z = as.numeric(Z),
    pval = as.numeric(pval),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    summary = summary,
    details = list(M1 = M1, M2 = M2, burn_in = burn_in, batch_size = batch_size)
  )
}
