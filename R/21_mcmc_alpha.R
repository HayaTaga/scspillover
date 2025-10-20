#' Horseshoe-Gibbs for synthetic weights alpha
#' @keywords internal
# hs_alpha_gibbs <- function(y, X, M = 2000, burn = 1000, verbose = TRUE) {
#   # y = X alpha + eps, alpha ~ HS(0); Makalic & Schmidt (2015) の補助変数化
#   stopifnot(is.numeric(y), is.matrix(X), length(y) == nrow(X))
#   T0 <- length(y)
#   N <- ncol(X)
#   XtX <- crossprod(X)
#   Xty <- crossprod(X, y)

#   alpha <- rep(0, N)
#   sigma2 <- var(y)
#   lambda2 <- rep(1, N)
#   tau2 <- 1
#   nu <- rep(1, N)
#   xi <- 1

#   draws <- matrix(NA_real_, nrow = M, ncol = N)
#   pb <- if (verbose) {
#     progress::progress_bar$new(total = M + burn, clear = FALSE)
#   } else {
#     NULL
#   }

#   for (iter in seq_len(M + burn)) {
#     if (verbose) {
#       pb$tick()
#     }

#     Dinv <- diag(1 / (lambda2 * tau2), N)
#     A <- XtX + Dinv
#     cholA <- chol(A)
#     mu <- backsolve(cholA, forwardsolve(t(cholA), Xty))
#     z <- rnorm(N)
#     alpha <- as.vector(mu + backsolve(cholA, z) * sqrt(sigma2))

#     res <- y - as.vector(X %*% alpha)
#     shape <- (T0 + N) / 2
#     rate <- (crossprod(res) + sum(alpha^2 / (lambda2 * tau2))) / 2
#     sigma2 <- 1 / rgamma(1, shape = shape, rate = rate)

#     for (j in 1:N) {
#       rate_l <- 1 / nu[j] + alpha[j]^2 / (2 * tau2 * sigma2)
#       lambda2[j] <- 1 / rgamma(1, shape = 1, rate = rate_l)
#       nu[j] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / lambda2[j])
#     }

#     rate_t <- 1 / xi + sum(alpha^2 / (2 * lambda2 * sigma2))
#     tau2 <- 1 / rgamma(1, shape = (N + 1) / 2, rate = rate_t)
#     xi <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / tau2)

#     if (iter > burn) draws[iter - burn, ] <- alpha
#   }
#   colnames(draws) <- paste0("alpha_", seq_len(N))
#   list(alpha = draws)
# }

hs_alpha_gibbs <- function(y, X, M = 2000, burn = 1000, verbose = TRUE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  draws <- hs_alpha_gibbs_cpp(y, X, M, burn, verbose)
  colnames(draws) <- paste0("alpha_", seq_len(ncol(X)))
  list(alpha = draws)
}
