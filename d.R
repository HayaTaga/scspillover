library(cmdstanr)
library(dplyr)
library(tidyr)
library(Matrix)

load("./data/california_smoking.rda")
california_smoking


source("R/01_utils.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/10_sc_spillover.R")
source("R/03_utils_plot.R")
source("R/40_geweke.R")
source("R/04_utils_diagnostics.R")
Rcpp::sourceCpp("./R/20_mcmc.cpp")
Rcpp::sourceCpp("./R/40_geweke.cpp")

panel_df <- california_smoking$panel
w_vec <- california_smoking$w
W_mat <- california_smoking$W
panel_df <- panel_df %>%
  mutate(
    treatment = ifelse((state == "California") & (year >= 1988), 1, 0)
  )

w <- as.matrix(w_vec[, 2])
W <- as.matrix(W_mat[, -1])

Yc_pre <- panel_df %>%
  filter(state_id != 0 & year <= 1989) %>%
  select(state, year, cigsale) %>%
  pivot_wider(
    names_from = state,
    values_from = cigsale,
    id_cols = year
  ) %>%
  select(-year)

Xc_pre <- panel_df %>%
  filter(state_id != 0 & year <= 1989) %>%
  select(state, year, retprice) %>%
  pivot_wider(
    names_from = state,
    values_from = retprice,
    id_cols = year
  ) %>%
  select(-year)

Yc_raw <- Yc_pre
T0 <- nrow(Yc_raw)
N <- ncol(Yc_raw)
K <- 1
p <- 1


X_list <- lapply(seq_len(T0), function(t) {
  X_array[t, , ] # これで N×K 行列になる
})

X_array <- array(0, dim = c(T0, N, K))
X_array[,, 1] <- as.matrix(Xc_pre)
X_list <- lapply(seq_len(T0), function(t) X_array[t, , ])

w_orig <- w
eigW <- eigen(W, only.values = TRUE)$values
maxabs <- max(abs(eigW))
bnd <- 0.95 / max(1.0, maxabs)

# ハイパーパラメータ
sigma_prior_scale <- 1.0
sigma_g_prior_scale <- 1.0
scale_global_alpha <- 1.0
scale_global_beta <- 1.0
nu_global <- 1.0
nu_local <- 1.0

stan_data <- list(
  T0 = T0,
  N = N,
  K = K,
  p = p,

  Yc_raw = Yc_raw,
  X = X_array,

  w_orig = w_orig,
  W = W_mat,

  bnd = bnd,

  scale_global_alpha = scale_global_alpha,
  scale_global_beta = scale_global_beta,
  nu_global = nu_global,
  nu_local = nu_local,
  sigma_prior_scale = sigma_prior_scale
)

library(cmdstanr)

# 1. コンパイル
mod <- cmdstan_model("R/20_mcmc.stan")

# 2. サンプリング
fit <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.9, # 必要に応じて上げる (0.9~0.99)
  max_treedepth = 12 # 必要に応じて上げる
)

print(fit, max_rows = 50)


library(haven)
df <- read_dta("./data/tabacco/synth_smoking.dta") %>%
  mutate(
    state_id = as.integer(state),
    state_name = as.character(haven::as_factor(state))
  ) %>%
  select(state_id, state_name, everything(), -state)

df
