## -----------------------------------------------------------------------------
rm(list = ls())
# library(scspill)
library(dplyr)
library(tidyr)
library(stringr)
library(knitr)
library(kableExtra)

run_mode <- tolower(Sys.getenv("SCSPILL_MODE", "full"))
if (!run_mode %in% c("full", "smoke")) {
  stop("SCSPILL_MODE must be either 'full' or 'smoke'.")
}
if (run_mode == "smoke") {
  m1_iter <- 5000L
  m2_iter <- 5000L
  burn_iter <- 2000L
} else {
  m1_iter <- 2000000L
  m2_iter <- 2000000L
  burn_iter <- 2000000L
}

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
if (!interactive()) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
}
seed_global <- 20251030L
set.seed(seed_global)
message(sprintf("[seed-init] script=geweke_main global_seed=%d mode=%s", seed_global, run_mode))

source("R/01_utils.R")
source("R/10_sc_spillover.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/30_simulation_.R")
source("R/40_geweke_latest.R")
source("R/41_robustness_check.R")
Rcpp::sourceCpp("src/20_mcmc.cpp")
Rcpp::sourceCpp("src/40_geweke_latest.cpp")


## -----------------------------------------------------------------------------
# --- dimensions ---
T0 <- 15
N <- 8
K <- 2 # X covariates
p <- 1 # latent factor

# --- W, w, alpha_hat_scaled ---
W_raw <- matrix(0, N, N)
for (i in 1:N) {
  if (i > 1) {
    W_raw[i, i - 1] <- 1
  }
  if (i < N) W_raw[i, i + 1] <- 1
}
w_raw <- rep(0, N)
w_raw[1] <- 1

W_use <- row_normalize(W_raw)
w_use <- as.numeric(w_raw)
wsum <- sum(w_use)
if (is.finite(wsum) && wsum > 1e-12) {
  w_use <- w_use / wsum
}

# fixed alpha_hat_scaled (in practice, scale BSCM alpha)
alpha_hat_scaled <- rnorm(N, 0, 0.4)

# --- X: T0 x N x K ---
Xc_pre <- array(0, dim = c(T0, N, K))
for (t in 1:T0) {
  for (i in 1:N) {
    for (k in 1:K) {
      Xc_pre[t, i, k] <- rnorm(1)
    }
  }
}

# --- true parameters (for internal consistency) ---
rho_true <- 0.25
sigma2_true <- 0.5^2
beta_true <- c(0.7, -0.5)
Eta_true <- matrix(rnorm(N * p, 0, 0.5), N, p)
Gamma_true <- matrix(0, p, T0)
phi_true <- 0.5
s2g_true <- 0.2
for (t in 1:T0) {
  if (t == 1) {
    Gamma_true[, t] <- rnorm(p, 0, sqrt(s2g_true / (1 - phi_true^2)))
  }
  if (t >= 2) {
    Gamma_true[, t] <- phi_true *
      Gamma_true[, t - 1] +
      rnorm(p, 0, sqrt(s2g_true))
  }
}

# --- prior ---
a0 <- 1.0
b0 <- 1.0 # IG(a0, b0) prior for sigma2
step_rho <- 0.05

# --- forward simulator in R (mirror of C++ forward kernel) ---
sim_forward_R <- function(
  T0,
  W_use,
  w_use,
  alpha_hat_scaled,
  rho,
  sigma2,
  Xc_pre,
  beta,
  Eta,
  Gamma
) {
  N <- nrow(W_use)
  K <- dim(Xc_pre)[3]
  A <- W_use + w_use %*% t(alpha_hat_scaled)
  I <- diag(N)
  Yc <- matrix(NA_real_, T0, N)
  for (t in 1:T0) {
    mu <- rep(0, N)
    if (K > 0) {
      Xt <- Xc_pre[t, , ] # N x K
      mu <- mu + as.vector(Xt %*% beta)
    }
    if (ncol(Eta) > 0) {
      mu <- mu + as.vector(Eta %*% Gamma[, t])
    }
    eps <- rnorm(N, 0, sqrt(sigma2))
    rhs <- mu + eps
    M <- I - rho * A
    Yc[t, ] <- as.vector(solve(M, rhs))
  }
  Yc
}

# pseudo observed treated pre-outcome (for g-function)
Y0_pre <- rnorm(T0)


## -----------------------------------------------------------------------------
# JDT parameters
M1 <- m1_iter
M2 <- m2_iter
burn_in <- burn_iter

out <- geweke_jdt_full(
  Y0_pre = Y0_pre,
  Yc_pre_like_dims = c(T0, N),
  W = W_use,
  w = w_use,
  alpha_hat_scaled = alpha_hat_scaled,
  Xc_pre = Xc_pre,
  p = p,
  M1 = M1,
  M2 = M2,
  burn_in = burn_in,
  a0 = a0,
  b0 = b0,
  step_rho = step_rho,
  g_fn = default_g_fn,
  verbose = TRUE
)

print(out$summary)

write.csv(out$summary, "output/tables/geweke_jdt_summary.csv", row.names = FALSE)


## -----------------------------------------------------------------------------
jdt_tbl <- out$summary

row_renamer <- c(
  rho = "$\\rho$",
  log_sigma2 = "$\\log(\\sigma^2)$",
  yc_mean = "$\\bar{y}_c$",
  log_yc_var = "$\\log(\\mathrm{Var}(y_c))$",
  spatial_quadratic = "$h_{\\text{spatial}}$",
  corr_y0_wyc = "$\\mathrm{Corr}(y_0, Wy)$",
  beta_mean = "$\\bar{\\beta}$",
  Eta_mean = "$\\bar{\\eta}$",
  Gamma_mean = "$\\bar{\\Gamma}$"
)

jdt_tbl$g <- unname(ifelse(
  jdt_tbl$g %in% names(row_renamer),
  row_renamer[jdt_tbl$g],
  jdt_tbl$g
))

colnames(jdt_tbl) <- c(
  "Statistic $g$",
  "Mean (iid)",
  "Mean (MCMC)",
  "SE (iid)",
  "SE (MCMC)",
  "$Z$",
  "$p$-value"
)

knitr::kable(jdt_tbl, caption = "Summary of Joint Distribution Test Results")

jdt_tex <- kableExtra::kable(
  jdt_tbl,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Summary of Joint Distribution Test Results"
)
cat(jdt_tex, file = "./output/tables/geweke_jdt_summary.tex")


## -----------------------------------------------------------------------------
thin_app <- 5L
out_app_thin <- geweke_jdt_full(
  Y0_pre = Y0_pre,
  Yc_pre_like_dims = c(T0, N),
  W = W_use,
  w = w_use,
  alpha_hat_scaled = alpha_hat_scaled,
  Xc_pre = Xc_pre,
  p = p,
  M1 = M1, # MC path unchanged
  M2 = M2 * thin_app, # run longer then thin inside g_fn if you implement that
  burn_in = burn_in,
  a0 = a0,
  b0 = b0,
  step_rho = 0.05,
  g_fn = default_g_fn,
  verbose = FALSE,
  rho_support = NULL
)

# (b) Custom g-function (example): centered squared norm of Y and rho^2
g_fn_custom <- function(theta, Y, W, w, Y0_pre) {
  yc <- as.numeric(Y)
  c(
    rho = theta$rho,
    rho_sq = theta$rho^2,
    log_sigma2 = log(theta$sigma2 + 1e-12),
    yc_mean = mean(yc),
    yc_center_q = mean((yc - mean(yc))^2)
  )
}

out_app_custom <- geweke_jdt_full(
  Y0_pre = Y0_pre,
  Yc_pre_like_dims = c(T0, N),
  W = W_use,
  w = w_use,
  alpha_hat_scaled = alpha_hat_scaled,
  Xc_pre = Xc_pre,
  p = p,
  M1 = M1,
  M2 = M2,
  burn_in = burn_in,
  a0 = a0,
  b0 = b0,
  step_rho = 0.05,
  g_fn = g_fn_custom, # <- custom statistics
  verbose = FALSE,
  rho_support = NULL
)

knitr::kable(
  out_app_custom$summary,
  caption = "Appendix JDT with custom g-function."
)
