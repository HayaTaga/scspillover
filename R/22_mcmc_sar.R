# # ===== 共通：ρのメトロポリス（LeSage型対数尤度） =====
# rho_loglik <- function(
#   rho,
#   Yc_pre,
#   Y0_pre,
#   Xc_pre,
#   beta,
#   W,
#   w,
#   Eta,
#   gamma,
#   sig2_e
# ) {
#   N <- ncol(Yc_pre)
#   T0 <- nrow(Yc_pre)
#   A <- diag(N) - rho * W
#   K <- if (is.null(Xc_pre)) 0 else dim(Xc_pre)[3]
#   logdetA <- as.numeric(determinant(A, logarithm = TRUE)$modulus)
#   ss <- 0
#   for (t in 1:T0) {
#     yc <- Yc_pre[t, ]
#     Xt <- Matrix(
#       matrix(as.numeric(Xc_pre[t, , , drop = FALSE]), ncol = K),
#       sparse = TRUE
#     )
#     xterm <- if (is.null(Xc_pre)) 0 else Xt %*% beta
#     u <- A %*%
#       yc -
#       rho * w * Y0_pre[t] -
#       as.matrix(xterm) -
#       as.matrix(Eta %*% gamma[, t])
#     ss <- ss + sum(u^2)
#   }
#   T0 * logdetA - (1 / (2 * sig2_e)) * ss
# }

# rho_mh_step <- function(
#   rho_cur,
#   Yc_pre,
#   Y0_pre,
#   Xc_pre,
#   beta,
#   W,
#   w,
#   Eta,
#   gamma,
#   sig2_e,
#   step = 0.05,
#   eigW_max = NULL
# ) {
#   if (is.null(eigW_max)) {
#     eigW_max <- max(abs(eigen(W, only.values = TRUE)$values))
#   }
#   bnd <- 0.95 / max(1, eigW_max)
#   prop <- rho_cur + rnorm(1, sd = step)
#   if (abs(prop) >= bnd) {
#     return(list(rho = rho_cur, accept = FALSE))
#   }
#   ll_cur <- rho_loglik(
#     rho_cur,
#     Yc_pre,
#     Y0_pre,
#     Xc_pre,
#     beta,
#     W,
#     w,
#     Eta,
#     gamma,
#     sig2_e
#   )
#   ll_prp <- rho_loglik(
#     prop,
#     Yc_pre,
#     Y0_pre,
#     Xc_pre,
#     beta,
#     W,
#     w,
#     Eta,
#     gamma,
#     sig2_e
#   )
#   if (log(runif(1)) < (ll_prp - ll_cur)) {
#     list(rho = prop, accept = TRUE)
#   } else {
#     list(rho = rho_cur, accept = FALSE)
#   }
# }

# # ===== Horseshoe（β用）の初期化と1ステップ更新 =====
# hs_beta_init <- function(P) {
#   list(
#     beta = rep(0, P),
#     sigma2 = 1,
#     lambda2 = rep(1, P),
#     tau2 = 1,
#     nu = rep(1, P),
#     xi = 1
#   )
# }

# hs_beta_update <- function(y, X, state) {
#   P <- ncol(X)
#   beta <- state$beta
#   sigma2 <- state$sigma2
#   lambda2 <- state$lambda2
#   tau2 <- state$tau2
#   nu <- state$nu
#   xi <- state$xi

#   Dinv <- diag(1 / (lambda2 * tau2), P)
#   XtX <- crossprod(X)
#   Xty <- crossprod(X, y)
#   A <- XtX + Dinv
#   cholA <- chol(A)
#   mu <- backsolve(cholA, forwardsolve(t(cholA), Xty))
#   z <- rnorm(P)
#   beta <- as.vector(mu + backsolve(cholA, z) * sqrt(sigma2))

#   res <- y - as.vector(X %*% beta)
#   shape <- (length(y) + P) / 2
#   rate <- (sum(res^2) + sum(beta^2 / (lambda2 * tau2))) / 2
#   sigma2 <- 1 / rgamma(1, shape = shape, rate = rate)

#   for (j in 1:P) {
#     rate_l <- 1 / nu[j] + (beta[j]^2) / (2 * tau2 * sigma2)
#     lambda2[j] <- 1 / rgamma(1, shape = 1, rate = rate_l)
#     nu[j] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / lambda2[j])
#   }

#   rate_t <- 1 / xi + sum(beta^2 / (2 * lambda2 * sigma2))
#   tau2 <- 1 / rgamma(1, shape = (P + 1) / 2, rate = rate_t)
#   xi <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / tau2)

#   list(
#     beta = beta,
#     sigma2 = sigma2,
#     lambda2 = lambda2,
#     tau2 = tau2,
#     nu = nu,
#     xi = xi
#   )
# }

# # ===== 潜在因子：FFBSとEta更新（簡潔版） =====
# ffbs_ar1 <- function(R_list, Eta, phi, sig2_g, sig2_e) {
#   T0 <- length(R_list)
#   p <- ncol(Eta)
#   N <- nrow(Eta)
#   H <- -Eta
#   Q <- sig2_g * diag(p)
#   Rm <- sig2_e * diag(N)
#   a <- matrix(0, p, T0)
#   P <- array(0, c(p, p, T0))
#   m <- matrix(0, p, T0)
#   C <- array(0, c(p, p, T0))
#   m0 <- rep(0, p)
#   C0 <- diag(1000, p)

#   for (t in 1:T0) {
#     if (t == 1) {
#       a[, t] <- phi * m0
#       P[,, t] <- phi * C0 %*% t(diag(p) * phi) + Q
#     } else {
#       a[, t] <- phi * m[, t - 1]
#       P[,, t] <- phi * C[,, t - 1] %*% t(diag(p) * phi) + Q
#     }
#     y_t <- R_list[[t]]
#     S <- H %*% P[,, t] %*% t(H) + Rm
#     K <- P[,, t] %*% t(H) %*% solve(S)
#     m[, t] <- a[, t] + K %*% (y_t - H %*% a[, t])
#     C[,, t] <- (diag(p) - K %*% H) %*% P[,, t]
#   }
#   gamma <- matrix(0, p, T0)
#   gamma[, T0] <- MASS::mvrnorm(1, mu = m[, T0], Sigma = C[,, T0])
#   for (t in (T0 - 1):1) {
#     if (t <= 0) {
#       break
#     }
#     B <- C[,, t] %*% t(diag(p) * phi) %*% solve(P[,, t + 1])
#     mean_t <- m[, t] + B %*% (gamma[, t + 1] - a[, t + 1])
#     cov_t <- C[,, t] - B %*% P[,, t + 1] %*% t(B)
#     gamma[, t] <- MASS::mvrnorm(1, mu = mean_t, Sigma = cov_t)
#   }
#   gamma
# }

# update_eta <- function(R_list, gamma, Sigma_eta, sig2_e) {
#   T0 <- length(R_list)
#   N <- length(R_list[[1]])
#   p <- nrow(gamma)
#   G <- t(gamma)
#   GG <- crossprod(G)
#   Eta <- matrix(0, N, p)
#   for (i in 1:N) {
#     r_i <- vapply(R_list, function(v) v[i], 0.0)
#     A <- GG / sig2_e + solve(Sigma_eta)
#     b <- crossprod(G, r_i) / sig2_e
#     cholA <- chol(A)
#     mu <- backsolve(cholA, forwardsolve(t(cholA), b))
#     z <- rnorm(p)
#     Eta[i, ] <- as.vector(mu + backsolve(cholA, z))
#   }
#   Eta
# }

# # ===== 本体：SAR + X + 因子 の MwG =====
# #' @keywords internal
# sar_gibbs_sampler <- function(
#   Y0_pre,
#   Yc_pre,
#   Xc_pre,
#   W,
#   w,
#   M = 2000,
#   burn = 1000,
#   verbose = TRUE,
#   p_factors = 1,
#   phi_gamma = 0.7,
#   sig2_g_init = 1.0,
#   sig2_e_init = 1.0,
#   step_rho = 0.05
# ) {
#   T0 <- nrow(Yc_pre)
#   N <- ncol(Yc_pre)
#   K <- if (is.null(Xc_pre)) 0 else dim(Xc_pre)[3]
#   p <- max(0, p_factors)

#   rho <- 0
#   sig2_e <- sig2_e_init
#   sig2_g <- sig2_g_init
#   Eta <- if (p > 0) matrix(0, N, p) else matrix(0, N, 0)
#   gamma <- if (p > 0) matrix(0, p, T0) else matrix(0, 0, T0)
#   beta_state <- if (K > 0) hs_beta_init(K) else NULL

#   Mtot <- M + burn
#   acc <- 0L
#   rho_draws <- numeric(M)
#   beta_draws <- if (K > 0) matrix(NA_real_, M, K) else NULL
#   eta_draws <- if (p > 0) array(NA_real_, c(N, p, M)) else NULL
#   gam_draws <- if (p > 0) array(NA_real_, c(p, T0, M)) else NULL
#   sig2e_draws <- numeric(M)

#   pb <- if (verbose) {
#     progress::progress_bar$new(total = Mtot, clear = FALSE)
#   } else {
#     NULL
#   }

#   for (it in 1:Mtot) {
#     if (verbose) {
#       pb$tick()
#     }

#     # R_t = (I - rho W) Yc_t - rho w Y0_t - X_t beta
#     A <- diag(N) - rho * W
#     R_list <- vector("list", T0)
#     for (t in 1:T0) {
#       yc <- as.matrix(Yc_pre[t, ])
#       Xt <- Matrix(
#         matrix(as.numeric(Xc_pre[t, , , drop = FALSE]), ncol = K),
#         sparse = TRUE
#       )
#       b <- Matrix(matrix(beta_state$beta, ncol = 1), sparse = TRUE)
#       xterm <- if (K > 0) Xt %*% b else rep(0, N)
#       R_list[[t]] <- as.vector(A %*% yc - rho * w * Y0_pre[t] - xterm)
#     }

#     if (p > 0) {
#       gamma <- ffbs_ar1(
#         R_list = R_list,
#         Eta = Eta,
#         phi = phi_gamma,
#         sig2_g = sig2_g,
#         sig2_e = sig2_e
#       )

#       R_list_eta <- vector("list", T0)
#       for (t in 1:T0) {
#         R_list_eta[[t]] <- R_list[[t]] + as.vector(Eta %*% gamma[, t])
#       }
#       Sigma_eta <- diag(10, p)
#       Eta <- update_eta(
#         R_list = R_list_eta,
#         gamma = gamma,
#         Sigma_eta = Sigma_eta,
#         sig2_e = sig2_e
#       )
#     }

#     if (K > 0) {
#       y_stack <- numeric(T0 * N)
#       X_stack <- matrix(0, T0 * N, K)
#       for (t in 1:T0) {
#         idx <- ((t - 1) * N + 1):(t * N)
#         y_stack[idx] <- if (p > 0) {
#           R_list[[t]] + as.vector(Eta %*% gamma[, t])
#         } else {
#           R_list[[t]]
#         }
#         X_stack[idx, ] <- Xc_pre[t, , ]
#       }
#       beta_state <- hs_beta_update(y = y_stack, X = X_stack, state = beta_state)
#     }

#     mh <- rho_mh_step(
#       rho_cur = rho,
#       Yc_pre = Yc_pre,
#       Y0_pre = Y0_pre,
#       Xc_pre = Xc_pre,
#       beta = if (K > 0) beta_state$beta else rep(0, K),
#       W = W,
#       w = w,
#       Eta = if (p > 0) Eta else matrix(0, N, 0),
#       gamma = if (p > 0) gamma else matrix(0, 0, T0),
#       sig2_e = sig2_e,
#       step = step_rho
#     )
#     rho <- mh$rho
#     if (mh$accept && it > burn) {
#       acc <- acc + 1L
#     }

#     # measurement variance
#     ss <- 0
#     for (t in 1:T0) {
#       u <- if (p > 0) {
#         R_list[[t]] - as.vector(Eta %*% gamma[, t])
#       } else {
#         R_list[[t]]
#       }
#       ss <- ss + sum(u^2)
#     }
#     shape <- (N * T0) / 2 + 0.1
#     rate <- ss / 2 + 0.1
#     sig2_e <- 1 / rgamma(1, shape = shape, rate = rate)

#     if (it > burn) {
#       k <- it - burn
#       rho_draws[k] <- rho
#       if (K > 0) {
#         beta_draws[k, ] <- beta_state$beta
#       }
#       if (p > 0) {
#         eta_draws[,, k] <- Eta
#         gam_draws[,, k] <- gamma
#       }
#       sig2e_draws[k] <- sig2_e
#     }
#   }

#   list(
#     rho = rho_draws,
#     beta = beta_draws,
#     Eta = eta_draws,
#     gamma = gam_draws,
#     sigma2_e = sig2e_draws,
#     acc_rate = acc / max(1, M)
#   )
# }

#' Full SAR MCMC (rho, beta, sigma2, Lambda, F) with optional covariates and factors
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
  step_rho = 0.05,
  c_beta = 10,
  c_lambda = 10,
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
    c_beta,
    c_lambda,
    a0,
    b0,
    verbose
  )

  list(
    rho = out$rho,
    beta = if (K > 0) out$beta else NULL,
    sigma2 = out$sigma2,
    Lambda = if (p > 0) out$Lambda else NULL, # N x p x M
    F = if (p > 0) out$F else NULL, # p x T0 x M
    acc_rate = out$acc_rate
  )
}
