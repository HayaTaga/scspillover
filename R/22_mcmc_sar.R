#' Full SAR MCMC (joint over rho & alpha, + optional covariates/factors)
#' @keywords internal
sar_gibbs_sampler <- function(
  Y0_pre,
  Yc_pre,
  Xc_pre,
  W,
  w,
  M = 2000,
  burn = 1000,
  verbose = TRUE,
  p_factors = 0,
  step_rho = 0.01,
  step_alpha = 0.01,
  a0 = 1,
  b0 = 1
) {
  Yc_pre <- as.matrix(Yc_pre)
  storage.mode(Yc_pre) <- "double"
  Y0_pre <- as.numeric(Y0_pre)
  W <- as.matrix(W)
  storage.mode(W) <- "double"
  w <- as.numeric(w)

  T0 <- nrow(Yc_pre)
  N <- ncol(Yc_pre)

  if (!is.null(Xc_pre)) {
    if (length(dim(Xc_pre)) != 3L) {
      stop("Xc_pre must be T0 x N x K array or NULL.")
    }
    if (dim(Xc_pre)[1] != T0 || dim(Xc_pre)[2] != N) {
      stop("Xc_pre dims mismatch.")
    }
    K <- dim(Xc_pre)[3]
    # C++ 側は (t*N + i)*K + k の一次元ベクトルとして受け取る
    Xvec <- as.numeric(aperm(Xc_pre, c(1, 2, 3)))
  } else {
    K <- 0L
    Xvec <- numeric()
  }
  p <- as.integer(max(0, p_factors))

  out <- sar_full_sampler_cpp(
    Y0_pre,
    Yc_pre,
    if (K > 0) Xvec else NULL,
    T0,
    N,
    K,
    p,
    w,
    W,
    M,
    burn,
    step_rho,
    step_alpha,
    a0,
    b0,
    verbose
  )

  list(
    rho = out$rho, # length M
    alpha = out$alpha, # M x N
    beta = if (K > 0) out$beta else NULL,
    sigma2 = out$sigma2,
    Lambda = if (p > 0) out$Lambda else NULL, # N x p x M
    F = if (p > 0) out$F else NULL, # p x T0 x M
    acc_rate = list(
      rho = out$acc_rho,
      alpha = out$acc_alpha
    )
  )
}
