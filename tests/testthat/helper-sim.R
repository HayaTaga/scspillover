rook_W <- function(r) {
  N <- r * r
  W <- matrix(0, N, N)
  id <- function(i, j) (i - 1) * r + j
  for (i in 1:r) for (j in 1:r) {
    k <- id(i, j)
    if (i > 1) W[k, id(i - 1, j)] <- 1
    if (i < r) W[k, id(i + 1, j)] <- 1
    if (j > 1) W[k, id(i, j - 1)] <- 1
    if (j < r) W[k, id(i, j + 1)] <- 1
  }
  W / rowSums(W)
}

gen_sim <- function(N = 36, T = 30, T0 = 20, rho = 0.3) {
  r <- sqrt(N)
  stopifnot(r == round(r))
  W <- rook_W(r)
  w <- c(rep(1, 4), rep(0, N - 4)); w <- w / sum(w)

  alpha <- rep(0, N)
  alpha[1] <- 0.5; alpha[2] <- -0.2; alpha[3:4] <- 0.4; alpha[5:10] <- 0.1 / 6

  Y0 <- numeric(T)
  Yc <- matrix(0, T, N)
  X <- matrix(rnorm(T * N), T, N)
  beta <- 1

  for (t in 1:T) {
    A <- diag(N) - rho * W
    e <- rnorm(N)
    rhs <- rho * w * Y0[max(1, t - 1)] + X[t, ] * beta + e
    yc0 <- solve(A, rhs)
    Yc[t, ] <- yc0
    Y0[t] <- sum(alpha * yc0) + rnorm(1, 0, 0.1)
  }

  for (t in (T0 + 1):T) Y0[t] <- Y0[t] + rnorm(1, 1, 1)

  times <- seq_len(T)
  panel <- data.frame(
    time = rep(times, each = N + 1),
    unit = c(rep("T0", T), paste0("C", rep(1:N, times = T))),
    y = c(Y0, as.vector(Yc))
  )

  list(panel = panel, T0 = T0, W = W, w = w)
}