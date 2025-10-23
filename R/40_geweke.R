# 40_geweke.R
# Geweke (2004) Joint Distribution Test (JDT) harness for your SAR+factors Bayesian model.
# This script orchestrates:
#  - Marginal–Conditional simulator:   theta ~ p(theta), y ~ p(y|theta)
#  - Successive–Conditional simulator: y ~ p(y|theta), theta' ~ q(theta|y)  (q = your MCMC, 1-sweep)
#  - Test functions g(theta,y), asymptotic variance estimates, and Z-statistics
#
# Usage (example):
#   Rcpp::sourceCpp("40_geweke.cpp")              # compiles prior + forward simulator
#   Rcpp::sourceCpp("your_mcmc_file.cpp")         # must export sar_full_sampler_cpp (and optionally hs_alpha_gibbs_cpp)
#   out <- geweke_jdt(Y0_pre, Yc_dims = c(T0,N), W, w, Xc = X_array, p = p, M1 = 20000, M2 = 20000)
#   print(out$summary)
#   # Inspect diagnostics:
#   str(out)
#
# The default g() is designed to catch typical implementation issues. You can pass your own g().

# ---------- Utilities ----------

build_initial_state <- function(prior, T0, N, K, p) {
  list(
    rho = prior$rho,
    sigma2 = prior$sigma2,
    beta = if (K > 0) as.numeric(prior$beta) else numeric(),
    Lambda = if (p > 0) prior$Lambda else matrix(0, N, 0),
    F = if (p > 0) prior$F else matrix(0, 0, T0),

    # Horseshoe（β）の補助
    sig2_b0 = if (K > 0) rep(1, K) else numeric(),
    nu_sig_b0 = if (K > 0) rep(1, K) else numeric(),
    tau2_b0 = 1,
    nu_tau_b0 = 1,

    # sigma^2 のハイパー
    nu_sigma2 = 1,

    # 因子側の補助（p=0 のときはダミー）
    phi_g = 0,
    s2_g = 1,
    nu_s2_g = 1,
    omega_k = if (p > 0) rep(1, p) else numeric(),
    nu_omega_k = if (p > 0) rep(1, p) else numeric(),
    s2_eta = 1,
    nu_s2_eta = 1
  )
}

# Flatten 3D X array (T0 x N x K) into ((t*N + i)*K + k) layout expected by C++
flatten_Xc <- function(Xc) {
  if (is.null(Xc)) {
    return(NULL)
  }
  stopifnot(length(dim(Xc)) == 3)
  T0 <- dim(Xc)[1]
  N <- dim(Xc)[2]
  K <- dim(Xc)[3]
  out <- numeric(T0 * N * K)
  for (t in seq_len(T0) - 1L) {
    for (i in seq_len(N) - 1L) {
      for (k in seq_len(K) - 1L) {
        idx <- (t * N + i) * K + k
        out[idx + 1L] <- Xc[t + 1L, i + 1L, k + 1L]
      }
    }
  }
  out
}

# Batch-means estimator for asymptotic variance (non-overlapping)
# Returns: list(tau2 = spectral variance estimate, var_mean = tau2 / M)
var_mcmc_batchmeans <- function(x, b = NULL) {
  x <- x[is.finite(x)]
  M <- length(x)
  if (is.null(b)) {
    b <- max(2L, floor(sqrt(M)))
  }
  a <- floor(M / b)
  if (a < 2L) {
    # Fallback: naive (anti-conservative) if chain too short
    return(list(tau2 = stats::var(x), var_mean = stats::var(x) / M))
  }
  # Use first a*b observations
  x2 <- x[seq_len(a * b)]
  bm <- colMeans(matrix(x2, nrow = b, ncol = a))
  tau2 <- b * stats::var(bm) # spectral variance estimator
  list(tau2 = as.numeric(tau2), var_mean = as.numeric(tau2) / (a * b))
}

# Default test function g(theta, y, ...): returns a *named* numeric vector
default_g_fn <- function(theta, Yc, Y0_pre, W, w) {
  out <- c()
  # Parameters
  out["rho"] <- unname(theta$rho)
  out["sigma2"] <- unname(theta$sigma2)
  if (!is.null(theta$beta) && length(theta$beta) > 0) {
    out["beta_L2_mean"] <- mean(theta$beta^2)
  }
  if (!is.null(theta$Lambda) && length(theta$Lambda) > 0) {
    out["Lambda_Frob2_over_np"] <- sum(theta$Lambda^2) /
      (nrow(theta$Lambda) * ncol(theta$Lambda))
  }
  # Data summaries
  yc_vec <- as.numeric(Yc)
  out["yc_mean"] <- mean(yc_vec)
  if (length(yc_vec) > 1) {
    out["yc_var"] <- stats::var(yc_vec)
  }
  # spatial quadratic form
  N <- ncol(Yc)
  T0 <- nrow(Yc)
  out["spatial_quadratic"] <- sum(diag(Yc %*% W %*% t(Yc))) / (N * T0)
  # corr(Y0, w' Yc_t)
  wyc <- as.numeric(Yc %*% as.numeric(w))
  if (stats::sd(Y0_pre) > 0 && stats::sd(wyc) > 0) {
    out["corr_y0_wyc"] <- stats::cor(as.numeric(Y0_pre), wyc)
  }
  out
}

robust_g_fn <- function(
  theta,
  Yc,
  Y0_pre,
  W,
  w,
  cap_quadratic = 1e6,
  eps = 1e-12
) {
  out <- c()

  # ---- parameter-side summaries (finite by construction) ----
  out["rho"] <- unname(theta$rho)

  s2 <- as.numeric(theta$sigma2)
  out["log_sigma2"] <- log(pmax(eps, s2))

  if (!is.null(theta$beta) && length(theta$beta) > 0) {
    out["log_beta_L2_mean"] <- log(pmax(eps, mean(theta$beta^2)))
  }

  if (!is.null(theta$Lambda) && length(theta$Lambda) > 0) {
    Lam <- as.matrix(theta$Lambda)
    np <- nrow(Lam) * ncol(Lam)
    if (np > 0) {
      out["log_Lambda_Frob2_over_np"] <- log(pmax(eps, sum(Lam^2) / np))
    }
  }

  # ---- data-side summaries (finite and robust) ----
  yc_vec <- as.numeric(Yc)
  out["yc_mean"] <- mean(yc_vec)

  v <- tryCatch(stats::var(yc_vec), error = function(e) NA_real_)
  out["log_yc_var"] <- log(pmax(eps, v))

  N <- ncol(Yc)
  T0 <- nrow(Yc)
  sq <- sum(diag(Yc %*% W %*% t(Yc))) / (N * T0) # per-time quadratic form averaged
  # cap to avoid domination by extremely rare huge draws
  sq <- max(min(sq, cap_quadratic), -cap_quadratic)
  out["spatial_quadratic_capped"] <- sq

  wyc <- as.numeric(Yc %*% as.numeric(w))
  if (stats::sd(Y0_pre) > 0 && stats::sd(wyc) > 0) {
    out["corr_y0_wyc"] <- stats::cor(as.numeric(Y0_pre), wyc)
  }

  out
}

# Extract theta from sar_full_sampler_cpp() output for M=1
extract_theta_from_posterior <- function(post, T0, N, K, p) {
  theta <- list(
    rho = as.numeric(post$rho[1L]),
    sigma2 = as.numeric(post$sigma2[1L]),
    beta = if (K > 0) as.numeric(post$beta[1L, , drop = TRUE]) else numeric(),
    Lambda = if (p > 0) {
      post$Lambda[,, 1L, drop = FALSE][,, 1L]
    } else {
      matrix(0, nrow = N, ncol = 0)
    },
    F = if (p > 0) {
      post$F[,, 1L, drop = FALSE][,, 1L]
    } else {
      matrix(0, nrow = 0, ncol = T0)
    }
  )
  theta
}

# Sanity checks for dimensions and inputs
.check_inputs <- function(Y0_pre, W, w, Xc, p) {
  stopifnot(is.numeric(Y0_pre), is.matrix(W), is.numeric(w))
  T0 <- length(Y0_pre)
  N <- nrow(W)
  stopifnot(ncol(W) == N, length(w) == N)
  if (!is.null(Xc)) {
    stopifnot(length(dim(Xc)) == 3L, dim(Xc)[1] == T0, dim(Xc)[2] == N)
  }
  if (!is.numeric(p) || p < 0) {
    stop("p must be a non-negative integer")
  }
  invisible(TRUE)
}

# Ensure required compiled symbols exist
.ensure_symbols <- function() {
  if (!exists("sample_prior_theta_cpp")) {
    stop(
      "C++ function 'sample_prior_theta_cpp' not found. Compile 40_geweke.cpp via Rcpp::sourceCpp."
    )
  }
  if (!exists("simulate_Yc_given_theta_cpp")) {
    stop(
      "C++ function 'simulate_Yc_given_theta_cpp' not found. Compile 40_geweke.cpp via Rcpp::sourceCpp."
    )
  }
  if (!exists("sar_full_sampler_cpp")) {
    stop(
      "Posterior sampler 'sar_full_sampler_cpp' not found. Compile your MCMC C++ (the file that exports sar_full_sampler_cpp)."
    )
  }
  invisible(TRUE)
}

# ---------- Main JDT runner ----------
# Arguments:
#  - Y0_pre: numeric length T0
#  - W: N x N matrix, w: length N
#  - Xc: optional T0 x N x K array (or NULL) for controls' regressors
#  - p: integer, number of latent factors
#  - M1, M2: lengths of the two simulators
#  - a0, b0: sigma^2 prior hyperparameters (Inv-Gamma)
#  - step_rho: proposal sd for rho in your MH kernel (passed through to sar_full_sampler_cpp)
#  - g_fn: function(theta, Yc, Y0_pre, W, w) -> named numeric vector (test statistics)
#  - batch_size: batch size for MCMC variance (default sqrt(M2))
#
# Returns: a list with fields
#  - summary: data.frame with g, mean_iid, mean_mcmc, se_iid, se_mcmc, Z, pval
#  - g_iid: matrix M1 x G
#  - g_mcmc: matrix M2 x G
#  - theta_examples: list(theta_iid = first theta from iid sim, theta_mcmc = last theta from mcmc sim)
#
geweke_jdt <- function(
  Y0_pre,
  W,
  w,
  Xc = NULL,
  p = 0L,
  M1 = 20000L,
  M2 = 20000L,
  a0 = 1.0,
  b0 = 1.0,
  step_rho = 0.01,
  g_fn = default_g_fn,
  batch_size = NULL,
  verbose = TRUE
) {
  .check_inputs(Y0_pre, W, w, Xc, p)
  .ensure_symbols()

  T0 <- length(Y0_pre)
  N <- nrow(W)
  K <- if (is.null(Xc)) 0L else dim(Xc)[3L]
  Xflat <- if (K > 0L) flatten_Xc(Xc) else NULL

  # ---------- Marginal–Conditional (iid) ----------
  if (verbose) {
    message("[JDT] Running marginal–conditional (iid) simulator...")
  }
  g_iid_list <- vector("list", M1)
  theta_first <- NULL
  for (m in seq_len(M1)) {
    th <- sample_prior_theta_cpp(T0, N, K, p, W, a0, b0)
    if (is.null(theta_first)) {
      theta_first <- th
    }
    Yc <- simulate_Yc_given_theta_cpp(
      Y0_pre,
      Xflat,
      T0,
      N,
      K,
      p,
      w,
      W,
      th$rho,
      th$sigma2,
      th$beta,
      th$Lambda,
      th$F
    )
    g_iid_list[[m]] <- g_fn(th, Yc, Y0_pre, W, w)
    if (verbose && (m %% max(1L, floor(M1 / 5))) == 0L) {
      message("  ... ", m, "/", M1)
    }
  }
  # align g vector names
  g_names <- names(g_iid_list[[1L]])
  G <- length(g_names)
  g_iid <- matrix(NA_real_, nrow = M1, ncol = G, dimnames = list(NULL, g_names))
  for (m in seq_len(M1)) {
    g_iid[m, ] <- g_iid_list[[m]][g_names]
  }

  # ---------- Successive–Conditional (MCMC kernel, 1-sweep each) ----------
  if (verbose) {
    message("[JDT] Running successive–conditional (true posterior kernel)...")
  }

  g_mcmc <- matrix(
    NA_real_,
    nrow = M2,
    ncol = G,
    dimnames = list(NULL, g_names)
  )

  # 初期の state は prior から（theta + 補助変数を一式）
  th0 <- sample_prior_theta_cpp(T0, N, K, p, W, a0, b0)
  state <- build_initial_state(th0, T0, N, K, p)

  for (m in seq_len(M2)) {
    # 1) y | theta（現在の state から生成）
    Yc <- simulate_Yc_given_theta_cpp(
      Y0_pre,
      Xflat,
      T0,
      N,
      K,
      p,
      w,
      W,
      state$rho,
      state$sigma2,
      if (K > 0) state$beta else numeric(),
      if (p > 0) state$Lambda else matrix(0, N, 0),
      if (p > 0) state$F else matrix(0, 0, T0)
    )

    # 2) theta' ~ 1ステップ posterior 遷移（初期値＝直前の state）
    state <- sar_full_one_step_cpp(
      Y0_pre,
      Yc,
      if (K > 0) Xflat else NULL,
      T0,
      N,
      K,
      p,
      w,
      W,
      state,
      step_rho = step_rho,
      a0 = a0,
      b0 = b0,
      verbose = FALSE
    )

    # 3) g(theta', y)
    th_cur <- list(
      rho = state$rho,
      sigma2 = state$sigma2,
      beta = if (K > 0) state$beta else numeric(),
      Lambda = if (p > 0) state$Lambda else matrix(0, N, 0),
      F = if (p > 0) state$F else matrix(0, 0, T0)
    )
    g_mcmc[m, ] <- g_fn(th_cur, Yc, Y0_pre, W, w)[g_names]

    if (verbose && (m %% max(1L, floor(M2 / 5))) == 0L) {
      message("  ... ", m, "/", M2)
    }
  }

  # ---------- Statistics ----------
  mean_iid <- colMeans(g_iid, na.rm = TRUE)
  mean_mcmc <- colMeans(g_mcmc, na.rm = TRUE)

  # iid variance of the mean (account for missing)
  n_iid <- colSums(is.finite(g_iid))
  se_iid <- sqrt(
    apply(g_iid, 2L, function(z) stats::var(z, na.rm = TRUE)) / pmax(1L, n_iid)
  )

  # MCMC variance of the mean (batch means; drop NA inside)
  se_mcmc <- vapply(
    seq_len(G),
    function(j) {
      var_mcmc_batchmeans(g_mcmc[, j], b = batch_size)$var_mean
    },
    numeric(1L)
  ) |>
    sqrt()

  Z <- (mean_iid - mean_mcmc) / sqrt(se_iid^2 + se_mcmc^2)
  pval <- 2 * stats::pnorm(-abs(Z))

  summary <- data.frame(
    g = g_names,
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
    g_iid = g_iid,
    g_mcmc = g_mcmc,
    theta_examples = list(theta_iid = theta_first, theta_mcmc = th),
    details = list(M1 = M1, M2 = M2, batch_size = batch_size)
  )
}

# ---------- Convenience wrapper with defaults ----------
# Constructs Xc if omitted (K=0), runs with a modest number of draws for smoke test.
geweke_jdt_quick <- function(
  Y0_pre,
  W,
  w,
  Xc = NULL,
  p = 0L,
  M1 = 2000L,
  M2 = 2000L,
  ...
) {
  geweke_jdt(Y0_pre, W, w, Xc, p, M1 = M1, M2 = M2, ...)
}
