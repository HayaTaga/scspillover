#' Metropolis sampler for SAR rho using pre-treatment periods
#' Likelihood: prod_t |I - rho W| * exp(-1/(2 s2) * ||(I - rho W)Yc_t - rho w Y0_t||^2)
#' with s2 profiled out. Uniform prior on rho in (-0.95, 0.95).
#' @keywords internal
rho_mh_sampler <- function(Y0, Yc, w, W, M = 2000, burn = 1000, step = 0.05, verbose = TRUE) {
  stopifnot(is.matrix(Yc), length(Y0) == nrow(Yc))
  N <- ncol(Yc); T0 <- nrow(Yc)
  eigW <- eigen(W, only.values = TRUE)$values
  # Stationarity region approx; keep inside (-1/max|eig|, 1/max|eig|)
  bnd <- 0.95/min(1, max(abs(eigW)))
  loglik <- function(rho) {
    if (abs(rho) >= bnd) return(-Inf)
    A <- diag(N) - rho * W
    logdet <- determinant(A, logarithm = TRUE)$modulus
    ss <- 0
    for (t in 1:T0) {
      yc <- Yc[t, ]
      u <- A %*% yc - rho * w * Y0[t]
      ss <- ss + sum(u^2)
    }
    # profile sigma^2 -> s2_hat = ss/(N*T0)
    ll <- as.numeric(T0 * logdet) - (N*T0/2) * log(ss/(N*T0)) - (N*T0)/2
    return(ll)
  }

  rho <- 0
  acc <- 0
  draws <- numeric(M)
  cur <- loglik(rho)
  pb <- if (verbose) progress::progress_bar$new(total = M + burn, clear = FALSE) else NULL
  for (iter in seq_len(M + burn)) {
    if (verbose) pb$tick()
    prop <- rho + rnorm(1, sd = step)
    lp <- loglik(prop)
    if (log(runif(1)) < (lp - cur)) {
      rho <- prop; cur <- lp
      if (iter > burn) acc <- acc + 1
    }
    if (iter > burn) draws[iter - burn] <- rho
  }
  list(rho = draws, acc_rate = acc/max(1, M))
}