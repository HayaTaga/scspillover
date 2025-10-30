# 41_robustness_check.R
# ---------------------------------------------------------------
# Prior sensitivity & Prior predictive utilities for the SC-spill model
# ---------------------------------------------------------------
# Notes:
# - Comments are in English (per the user's request).
# - This file assumes the C++ kernel `scspill_one_step_cpp()` is available
#   with the signature that accepts `rho_lo` and `rho_hi`.
# - It also assumes `normalize_joint_wW()` and `compute_bnd()` exist in R.

# ===============================================================
# Posterior MCMC (data fixed) for prior sensitivity analysis
# ===============================================================
run_mcmc_for_posterior <- function(
  Yc_obs,               # T0 x N matrix (observed pre-period control outcomes)
  W_raw, w_raw,
  alpha_hat_scaled,     # length-N vector (BSCM alpha aligned to column scaling)
  Xc_pre = NULL,        # T0 x N x K array or NULL
  p = 0L,               # number of factors (0 is recommended initially)
  a0 = 1.0, b0 = 1.0,   # sigma^2 ~ InvGamma(a0, b0)
  rho_support = NULL,   # c(lo, hi); if NULL, a spectral bound is used
  step_rho = 0.05,
  M_burn = 5000L, M_keep = 20000L,
  thin = 1L,            # thinning interval (1 = no thinning)
  seed = 123
){
  # --- basic checks
  stopifnot(is.matrix(Yc_obs))
  set.seed(seed)

  T0 <- nrow(Yc_obs); N <- ncol(Yc_obs)
  K  <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  # --- normalize (w | W) row-wise, consistent with the model
  nm <- normalize_joint_wW(W_raw, w_raw)
  W  <- nm$W_new
  w  <- nm$w_new

  # --- ensure X is a 3D array (T0 x N x K), use empty array if K=0
  X_use <- if (is.null(Xc_pre)) {
    array(0, c(T0, N, 0L))
  } else {
    stopifnot(dim(Xc_pre)[1] == T0, dim(Xc_pre)[2] == N)
    Xc_pre
  }

  # --- support for rho
  if (is.null(rho_support)) {
    bnd <- compute_bnd(W, c_stability = 0.95)
    rho_lo <- -bnd; rho_hi <-  bnd
  } else {
    stopifnot(length(rho_support) == 2L, rho_support[1] < rho_support[2])
    rho_lo <- as.numeric(rho_support[1]); rho_hi <- as.numeric(rho_support[2])
  }

  # --- initialize state consistent with priors
  rinvgamma1 <- function(a, b) 1/rgamma(1, shape = a, rate = b)
  state <- list(
    rho    = as.numeric(runif(1, rho_lo, rho_hi)),
    sigma2 = as.numeric(rinvgamma1(a0, b0)),
    beta   = if (K > 0) rnorm(K, 0, 1) else numeric(0),
    Eta    = if (p > 0) matrix(0, N, p) else matrix(0, N, 0),
    Gamma  = if (p > 0) matrix(0, p, T0) else matrix(0, 0, T0)
  )

  # --- burn-in: apply a single-step kernel with observed data fixed
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

  # --- keep: thinning-aware collection
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

  # --- summarize posterior (extend as needed for beta, Eta, Gamma)
  rho_vec <- vapply(draws, function(s) s$rho,    numeric(1))
  s2_vec  <- vapply(draws, function(s) s$sigma2, numeric(1))

  out_df <- data.frame(
    param = c("rho","sigma2"),
    mean  = c(mean(rho_vec), mean(s2_vec)),
    sd    = c(sd(rho_vec),   sd(s2_vec)),
    q025  = c(unname(quantile(rho_vec, 0.025)), unname(quantile(s2_vec, 0.025))),
    q975  = c(unname(quantile(rho_vec, 0.975)), unname(quantile(s2_vec, 0.975)))
  )

  # beta summary if present (vector of length K)
  if (K > 0) {
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

  list(
    theta = out_df,
    raw   = list(rho = rho_vec, sigma2 = s2_vec),
    meta  = list(M_burn = M_burn, M_keep = M_keep, thin = thin,
                 rho_lo = rho_lo, rho_hi = rho_hi, step_rho = step_rho)
  )
}

# ===============================================================
# Prior sensitivity wrapper
# ===============================================================
prior_sensitivity <- function(
  Yc_obs, W_raw, w_raw,
  alpha_hat_scaled,
  Xc_pre = NULL, p = 0L,
  grid,                       # data.frame with columns: a0, b0, rho_lo, rho_hi, step_rho
  M_burn = 5000L, M_keep = 20000L,
  thin = 1L
){
  # Validate grid
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

# ===============================================================
# Pure R forward simulator (used in prior predictive checks)
# ===============================================================
simulate_Yc_forward_R <- function(
  T0, W_use, w_use, alpha_hat_scaled,
  rho, sigma2,
  Xc_pre, beta, Eta, Gamma
){
  # Shapes
  N <- nrow(W_use)
  K <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  # Build A = W + w * alpha^T  (use tcrossprod for clarity)
  A <- W_use + tcrossprod(w_use, as.numeric(alpha_hat_scaled))
  I <- diag(N)

  Yc <- matrix(NA_real_, T0, N)
  sd_eps <- sqrt(max(sigma2, 1e-12))

  for (t in seq_len(T0)) {
    mu <- rep(0, N)

    # Add X_t beta if present
    if (K > 0) {
      arr <- Xc_pre[t, , , drop = FALSE]       # 1 x N x K
      Nloc <- dim(arr)[2]; Kloc <- dim(arr)[3]
      Xt <- matrix(arr, nrow = Nloc, ncol = Kloc)  # N x K
      mu <- mu + as.numeric(Xt %*% beta)
    }

    # Add factor component if present
    if (ncol(Eta) > 0) mu <- mu + as.numeric(Eta %*% Gamma[, t])

    # Draw noise and solve (I - rho A) Y_t = mu + eps
    eps <- rnorm(N, 0, sd_eps)
    rhs <- mu + eps
    M <- I - rho * A

    Yt <- tryCatch(
      solve(M, rhs),
      error = function(e) {
        # Stable fallback: solve (M'M) Y = M' rhs
        AtA <- crossprod(M)
        diag(AtA) <- diag(AtA) + 1e-10
        solve(AtA, crossprod(M, rhs))
      }
    )
    Yc[t, ] <- Yt
  }
  Yc
}

# ===============================================================
# Full prior sampler for theta = {rho, sigma2, beta, Eta, Gamma}
# ===============================================================
prior_sampler_theta_full <- function(
  T0, N, K, p, a0, b0,
  W_use, w_use, alpha_hat_scaled,
  rho_support = NULL
){
  # rho ~ Unif(rho_lo, rho_hi) or spectral bound by default
  if (is.null(rho_support)) {
    bnd <- compute_bnd(W_use, c_stability = 0.95)
    rho <- runif(1, -bnd, bnd)
  } else {
    rho <- runif(1, rho_support[1], rho_support[2])
  }

  sigma2 <- 1/rgamma(1, shape = a0, rate = b0)
  beta   <- if (K > 0) rnorm(K, 0, 1) else numeric(0)
  Eta    <- if (p > 0) matrix(rnorm(N * p, 0, 1), N, p) else matrix(0, N, 0)
  Gamma  <- if (p > 0) matrix(rnorm(p * T0, 0, 1), p, T0) else matrix(0, 0, T0)

  list(rho = rho, sigma2 = sigma2, beta = beta, Eta = Eta, Gamma = Gamma)
}

# ===============================================================
# Prior predictive summary statistics
# ===============================================================
ppc_stats <- function(Yc, Y0_pre, W_use, w_use){
  # Compute a small set of informative discrepancy measures.
  yc <- as.numeric(Yc)
  N <- ncol(Yc); T0 <- nrow(Yc)

  wyc <- as.numeric(Yc %*% as.numeric(w_use))
  spatial_q <- sum(diag(Yc %*% W_use %*% t(Yc))) / (N * T0)

  ac1 <- tryCatch({
    # Mean-centered across time to compute average AR(1) by unit
    Yd <- scale(t(Yc), center = TRUE, scale = FALSE)  # N x T0
    num <- rowSums(Yd[, -1, drop = FALSE] * Yd[, -ncol(Yd), drop = FALSE])
    den <- rowSums(Yd * Yd)
    mean(num / den, na.rm = TRUE)
  }, error = function(e) NA_real_)

  c(
    yc_mean = mean(yc),
    log_yc_var = log(var(yc) + 1e-12),
    spatial_quadratic = spatial_q,
    corr_y0_wyc = if (sd(Y0_pre) > 0 && sd(wyc) > 0) cor(Y0_pre, wyc) else NA_real_,
    ac1 = ac1
  )
}

# ===============================================================
# Prior predictive main routine
# ===============================================================
prior_predictive <- function(
  Y0_pre, Yc_obs = NULL,       # if provided, observed stats are overlaid
  W_raw, w_raw, alpha_hat_scaled,
  Xc_pre = NULL, p = 0L,
  a0 = 3, b0 = 1, rho_support = NULL,
  R = 2000L, seed = 123
){
  set.seed(seed)

  T0 <- length(Y0_pre)
  N  <- nrow(W_raw)
  K  <- if (is.null(Xc_pre)) 0L else dim(Xc_pre)[3]

  # Normalize (w | W)
  nm <- normalize_joint_wW(W_raw, w_raw)
  W  <- nm$W_new
  w  <- nm$w_new

  # Observed statistics (optional)
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
      if (is.null(Xc_pre)) array(0, c(T0, N, 0)) else Xc_pre,
      th$beta, th$Eta, th$Gamma
    )
    stat_mat[r, ] <- ppc_stats(Yc_sim, Y0_pre, W, w)
  }

  list(
    stat = as.data.frame(stat_mat),
    observed = obs_stat
  )
}

# ===============================================================
# Optional: simple plotting helper for prior predictive histograms
# (Overlays a red vertical line for observed statistics if provided)
# ===============================================================
prior_predictive_plot <- function(pp_out, main_prefix = "Prior predictive"){
  stopifnot(is.list(pp_out), "stat" %in% names(pp_out))
  S <- pp_out$stat
  obs <- pp_out$observed

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(2, 3), mar = c(4,4,3,1))
  cols <- names(S)

  for (nm in cols) {
    h <- hist(S[[nm]], breaks = "FD", col = "grey85",
              border = "white", main = sprintf("%s: %s", main_prefix, nm),
              xlab = nm)
    if (!is.null(obs) && is.finite(obs[nm])) {
      abline(v = obs[nm], col = "red", lwd = 2)
    }
  }
  invisible(NULL)
}