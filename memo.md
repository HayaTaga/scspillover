# scspill — Synthetic Control with Spillovers (Bayesian)

このキャンバスには、論文「Identification and Inference for Synthetic Control Methods with Spillover Effects」（以下、論文）に基づく再現性パッケージの初期実装一式（Rパッケージ雛形＋コア実装＋ビネット＋テスト）を収めます。**R CMD build / check** を前提に、`roxygen2` によるドキュメント生成、`testthat` による自動テスト、`pkgdown` によるサイト化を想定した構成。

---

## 0. ディレクトリ構成（概要）

```
scspill/
├─ DESCRIPTION
├─ NAMESPACE                      # roxygen2 により自動生成
├─ R/
│   ├─ sc_spillover.R             # メインAPI（推定）
│   ├─ hs_alpha.R                 # horseshoeによるα推定（Gibbs）
│   ├─ rho_sar.R                  # ρの事後サンプリング（Metropolis）
│   ├─ effects.R                  # 式(5)(6)の計算、要約・可視化S3
│   ├─ utils_data_prep.R          # 前処理、行列生成
│   └─ zzz.R                      # パッケージ起動時メッセージ
├─ data/
│   ├─ cali_tobacco.rda          # カリフォルニア事例データ
│   └─ sudan_split.rda           # スーダン分裂データ
├─ man/                           # roxygen2が生成
├─ vignettes/
│   ├─ vignette-california.Rmd
│   ├─ vignette-sudan.Rmd
│   └─ vignette-simulation.Rmd
├─ tests/
│   └─ testthat/
│       ├─ test-alpha.R
│       ├─ test-rho.R
│       ├─ test-effects.R
│       └─ helper-sim.R
├─ README.md
└─ .Rbuildignore
```

---

## 1. DESCRIPTION

```r
Package: scspill
Title: Synthetic Control with Spillover Effects (Bayesian)
Version: 0.1.0
Authors@R: c(
    person(given = "Shosei", family = "Sakaguchi", role = c("aut", "cph")),
    person(given = "Hayato", family = "Tagawa", role = c("aut")),
    person(given = "Your Name", family = "(Maintainer)", email = "you@example.com", role = c("cre"))
)
Description: Implements a Bayesian synthetic control method that allows for
    spillover effects via a spatial autoregressive (SAR) structure. It estimates
    synthetic weights with horseshoe priors and samples the SAR parameter rho,
    then identifies treatment and spillover effects using analytic formulas in the
    paper. Includes vignettes reproducing California tobacco and Sudan split studies.
License: MIT + file LICENSE
Depends: R (>= 4.2.0)
Imports:
    stats,
    methods,
    Matrix,
    progress,
    ggplot2
Suggests:
    knitr,
    rmarkdown,
    testthat (>= 3.1.0),
    covr,
    pkgdown
Encoding: UTF-8
LazyData: true
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.1
VignetteBuilder: knitr
Config/testthat/edition: 3
```

---

## 2. NAMESPACE（roxygen2生成想定、参考）

```r
export(sc_spillover)
export(plot.scspill)
export(summary.scspill)
import(Matrix)
importFrom(stats, rnorm, runif, sd)
importFrom(ggplot2, ggplot, aes, geom_line, geom_ribbon, theme_minimal, labs)
S3method(plot, scspill)
S3method(summary, scspill)
```

---

## 3. R/sc\_spillover.R（メインAPI）

```r
#' Synthetic Control with Spillovers (Bayesian)
#'
#' @description
#' Implements the estimator proposed in the paper by combining
#' (i) horseshoe-Bayesian estimation of synthetic weights (alpha),
#' (ii) Metropolis sampling of spatial autoregressive parameter (rho), and
#' (iii) closed-form identification of treatment and spillover effects via Eqs. (5)(6).
#'
#' @param data data.frame with columns: unit, time, y (outcome). Optional: covariates X_*
#' @param treated_unit single value matching `unit` for the treated unit (i=0 in paper).
#' @param T0 integer; number of pre-treatment periods.
#' @param w numeric vector length N (weights from each control unit to treated unit).
#' @param W numeric NxN matrix (row-normalized spatial weight among control units).
#' @param M integer; number of posterior draws.
#' @param burn integer; burn-in iterations for Gibbs / MH.
#' @param seed integer; RNG seed.
#' @param verbose logical; progress display.
#'
#' @returns An object of class `scspill` with elements:
#' * alpha_draws (M x N), rho_draws (M), alpha_hat, rho_hat
#' * effects: list with posterior means and credible intervals for
#'   - treat: T-T0 length vector for treated unit effects (xi_0t)
#'   - spill: list of N control-unit time series (xi_it)
#' * inputs: list of matrices used (Y0, Yc_pre, Yc_post, etc.)
#' @export
sc_spillover <- function(data, treated_unit, T0, w, W, M = 2000, burn = 1000,
                         seed = 123, verbose = TRUE) {
  stopifnot(is.data.frame(data))
  set.seed(seed)

  prep <- scspill_prep(data, treated_unit = treated_unit, T0 = T0)
  Y0_pre  <- prep$Y0_pre            # length T0
  Yc_pre  <- prep$Yc_pre            # T0 x N
  Y0_post <- prep$Y0_post           # length T1
  Yc_post <- prep$Yc_post           # T1 x N
  N <- ncol(Yc_pre); T1 <- nrow(Yc_post)

  # (i) Horseshoe for alpha using pre-treatment: minimize (2)
  hs <- hs_alpha_gibbs(y = Y0_pre, X = Yc_pre, M = M, burn = burn, verbose = verbose)
  alpha_draws <- hs$alpha
  alpha_hat <- colMeans(alpha_draws)

  # (ii) Sample rho using pre-treatment SAR likelihood (no X for identification need)
  rho_out <- rho_mh_sampler(Y0 = Y0_pre, Yc = Yc_pre, w = w, W = W,
                            M = M, burn = burn, verbose = verbose)
  rho_draws <- rho_out$rho
  rho_hat <- mean(rho_draws)

  # (iii) Effects via identification (5)(6) using post-treatment observations only
  eff <- posterior_effects(Y0_post = Y0_post, Yc_post = Yc_post,
                           alpha_draws = alpha_draws, rho_draws = rho_draws,
                           w = w, W = W)

  structure(list(alpha_draws = alpha_draws,
                 rho_draws = rho_draws,
                 alpha_hat = alpha_hat,
                 rho_hat = rho_hat,
                 effects = eff,
                 inputs = prep),
            class = "scspill")
}
```

---

## 4. R/hs\_alpha.R（Horseshoeによる α 推定）

```r
#' Horseshoe-Gibbs for synthetic weights alpha
#' @keywords internal
hs_alpha_gibbs <- function(y, X, M = 2000, burn = 1000, verbose = TRUE) {
  # Model: y = X alpha + eps, alpha ~ HS(0), sigma^2 unknown
  # Makalic & Schmidt (2015) parameter augmentation
  T0 <- length(y); N <- ncol(X)
  XtX <- crossprod(X); Xty <- crossprod(X, y)

  # priors/initial
  alpha <- rep(0, N)
  sigma2 <- var(y)
  lambda2 <- rep(1, N)
  tau2 <- 1
  nu <- rep(1, N)   # for lambda2
  xi <- 1           # for tau2

  draws <- matrix(NA_real_, nrow = M, ncol = N)

  pb <- if (verbose) progress::progress_bar$new(total = M + burn, clear = FALSE) else NULL
  for (iter in seq_len(M + burn)) {
    if (verbose) pb$tick()

    # alpha | rest  ~ N( (XtX + D^-1)^-1 X'y, ... ) where D = diag(lambda2 * tau2)
    Dinv <- diag(1/(lambda2 * tau2), N)
    A <- XtX + Dinv
    cholA <- chol(A)
    mu <- backsolve(cholA, forwardsolve(t(cholA), Xty))
    z <- rnorm(N)
    alpha <- mu + backsolve(cholA, z) * sqrt(sigma2)

    # sigma2 | rest  ~ IG
    res <- y - as.vector(X %*% alpha)
    shape <- (T0 + N)/2
    rate <- (crossprod(res) + sum(alpha^2/(lambda2 * tau2)))/2
    sigma2 <- 1/rgamma(1, shape = shape, rate = rate)

    # lambda2 | rest  (Half-Cauchy via IG augmentation)
    for (j in 1:N) {
      rate_l <- 1/nu[j] + alpha[j]^2/(2 * tau2 * sigma2)
      lambda2[j] <- 1/rgamma(1, shape = 1, rate = rate_l)
      nu[j] <- 1/rgamma(1, shape = 1, rate = 1 + 1/lambda2[j])
    }

    # tau2 | rest
    rate_t <- 1/xi + sum(alpha^2/(2 * lambda2 * sigma2))
    tau2 <- 1/rgamma(1, shape = (N + 1)/2, rate = rate_t)
    xi <- 1/rgamma(1, shape = 1, rate = 1 + 1/tau2)

    if (iter > burn) draws[iter - burn, ] <- alpha
  }
  colnames(draws) <- paste0("alpha_", seq_len(N))
  list(alpha = draws)
}
```

---

## 5. R/rho\_sar.R（ρ のサンプリング：Metropolis）

```r
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
```

---

## 6. R/effects.R（式(5)(6)の実装、要約・可視化）

```r
#' Posterior effects via identification formulas (5)(6)
#' @keywords internal
posterior_effects <- function(Y0_post, Yc_post, alpha_draws, rho_draws, w, W, cred = 0.95) {
  T1 <- length(Y0_post); N <- ncol(Yc_post); M <- nrow(alpha_draws)
  Aeff <- array(NA_real_, dim = c(T1, M))          # treatment effects per draw
  Spill <- array(NA_real_, dim = c(T1, N, M))      # spillovers per draw

  IN <- diag(N)
  for (m in 1:M) {
    a <- alpha_draws[m, ]; r <- rho_draws[m]
    Ainv <- solve(IN - r * (outer(w, a) + W)) # (IN - r w a' - r W)^{-1}
    B <- (IN - r * W)                         # (IN - r W)
    for (t in 1:T1) {
      yc <- Yc_post[t, ]; y0 <- Y0_post[t]
      tmp <- Ainv %*% (B %*% yc - r * w * y0)
      # (5) treatment
      Aeff[t, m] <- y0 - as.numeric(crossprod(a, tmp))
      # (6) spillovers
      Spill[t, , m] <- yc - as.vector(tmp)
    }
  }
  treat_mean <- rowMeans(Aeff)
  treat_q <- apply(Aeff, 1, quantile, probs = c((1-cred)/2, 1 - (1-cred)/2))

  spill_mean <- apply(Spill, c(1,2), mean)
  spill_lo <- apply(Spill, c(1,2), quantile, probs = (1-cred)/2)
  spill_hi <- apply(Spill, c(1,2), quantile, probs = 1 - (1-cred)/2)

  list(treat = list(mean = treat_mean, lo = treat_q[1,], hi = treat_q[2,]),
       spill = list(mean = spill_mean, lo = spill_lo, hi = spill_hi))
}

#' @export
summary.scspill <- function(object, ...) {
  out <- list(rho_mean = mean(object$rho_draws),
              rho_ci = quantile(object$rho_draws, c(0.025, 0.975)),
              alpha_nonzero = which(abs(object$alpha_hat) > 1e-4),
              treat_avg = mean(object$effects$treat$mean))
  class(out) <- "summary.scspill"
  out
}

#' @export
print.summary.scspill <- function(x, ...) {
  cat("scspill summary\n")
  cat(sprintf("rho (post. mean) = %.3f, 95%% CI [%.3f, %.3f]\n",
              mean(x$rho_mean), 
```
