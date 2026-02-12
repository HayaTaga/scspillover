## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 4,
  message = FALSE,
  warning = FALSE
)


## -----------------------------------------------------------------------------
rm(list = ls())
# library(scspill)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(knitr)
library(kableExtra)
library(ggpattern)

run_mode <- tolower(Sys.getenv("SCSPILL_MODE", "full"))
if (!run_mode %in% c("full", "smoke")) {
  stop("SCSPILL_MODE must be either 'full' or 'smoke'.")
}
if (run_mode == "smoke") {
  mcmc_iter <- 8000L
  burn_iter <- 2000L
  sens_burn <- 5000L
  sens_keep <- 10000L
  ppc_rep <- 5000L
} else {
  mcmc_iter <- 1000000L
  burn_iter <- 500000L
  sens_burn <- 500000L
  sens_keep <- 1000000L
  ppc_rep <- 100000L
}

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
if (!interactive()) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
}
seed_global <- 20251022L
set.seed(seed_global)
message(sprintf(
  "[seed-init] script=sudan_main global_seed=%d mode=%s",
  seed_global,
  run_mode
))

source("R/01_utils.R")
source("R/02_utils_data_prep.R")
source("R/03_utils_plot.R")
source("R/04_utils_diagnostics.R")
source("R/10_sc_spillover.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/40_geweke_latest.R")
source("R/41_robustness_check.R")
source("R/50_plots_paper.R")
Rcpp::sourceCpp("src/20_mcmc.cpp")
Rcpp::sourceCpp("src/40_geweke_latest.cpp")


load("./data/sudan_secession.rda")

panel_df <- if ("panel" %in% names(sudan_secession)) {
  sudan_secession[["panel"]]
} else if ("panel_df" %in% names(sudan_secession)) {
  sudan_secession[["panel_df"]]
} else {
  stop("sudan_secession must contain 'panel' or 'panel_df'.")
}
w_vec <- sudan_secession$w
W_mat <- sudan_secession$W
panel_df <- panel_df %>%
  mutate(
    treatment = ifelse((country == "Sudan") & (year >= 2011), 1, 0)
  )


## -----------------------------------------------------------------------------
w <- as.matrix(w_vec[, 2])
W <- as.matrix(W_mat[, -1])

fit <- sc_spillover(
  data = panel_df,
  treated_unit = "Sudan",
  w = w,
  W = W,
  treatment_dummy = "treatment",
  y = colnames(panel_df)[4],
  X = colnames(panel_df)[c(5:10)],
  p_factors = 1,
  M = mcmc_iter,
  burn = burn_iter,
  seed = 20251022,
  step_rho = 0.02,
  unit_col = "country",
  time_col = "year"
)


## ----fig3b--------------------------------------------------------------------
plot(fit, time_col = "year")
plot(fit, type = "spill_top", top_n = 8, time_col = "year")
plot(fit, type = "rho")

p <- diagnostics.scspill(
  fit,
  which_alpha = NULL,
  top_n_alpha = 6,
  which_beta = NULL,
  top_n_beta = 6
)
print(p)

ggsave(
  filename = "output/figures/fig_sudan_diag.png",
  plot = p,
  width = 12.0,
  height = 6.0,
  dpi = 300,
)

diag_tab <- attr(p, "summary")
diag_tab[order(diag_tab$ess, decreasing = TRUE), ]

write.csv(diag_tab, "output/tables/sudan_diag_table.csv", row.names = FALSE)

digits_vec <- rep(3L, ncol(diag_tab))
names(digits_vec) <- names(diag_tab)
# digits_vec[c("n")] <- 0L
digits_vec[c(
  "mean",
  "sd",
  "q025",
  "q50",
  "q975",
  "ess",
  "rhat\\_split"
)] <- 4L
digits_vec <- digits_vec[colnames(diag_tab)]

col_labels <- c(
  parameter = "Parameter",
  mean = "Mean",
  sd = "SD",
  q025 = "$2.5\\%$",
  q50 = "$50\\%$",
  q975 = "$97.5\\%$",
  ess = "ESS",
  rhat_split = "$\\hat{R}$"
)
col_labels <- col_labels[colnames(diag_tab)]

latex_code <- kableExtra::kable(
  diag_tab,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  digits = as.integer(digits_vec),
  col.names = unname(col_labels),
  caption = "MCMC Posterior Summary and Diagnostics (Sudan Split)",
  label = "tab:sudan_diag"
) |>
  kableExtra::kable_styling(latex_options = c("hold_position", "scale_down"))

cat(latex_code, file = "output/tables/sudan_diag_table.tex")


## -----------------------------------------------------------------------------
print(round(fit$rho_hat, digits = 3))
print(round(quantile(fit$rho_draws, probs = c(0.025, 0.975)), digits = 3))

p_a <- plot_panel_outcomes(
  fit = fit,
  data = panel_df,
  treated_unit = "Sudan",
  treatment_dummy = "treatment",
  y = colnames(panel_df)[4],
  unit_col = "country",
  time_col = "year",
  title = "",
  x_title = "Year",
  y_title = colnames(panel_df)[4],
  year_every = 2, # 2年ごと表示
  force_perfect_pre = TRUE,
  skip_year = 2011,
  save_path = "output/figures/fig_sudan_panelA.pdf"
)

p_a

annote_main <- data.frame(
  time = 2010, # 矢印の付け根のX座標（年）
  y = -100, # 矢印の付け根のY座標
  label = "Secession",
  hjust = 1, # テキストを矢印の左側に配置
  vjust = 0.5
)

p_b <- plot_panel_effect(
  fit = fit,
  data = panel_df, # ★ SCM計算のために必須
  treated_unit = "Sudan",
  treatment_dummy = "treatment",
  y = colnames(panel_df)[4],
  unit_col = "country",
  time_col = "year",
  vline_at = 2010, # 論文の図の垂直線
  annotations_df = annote_main,
  skip_year = 2011
)

p_b

ggsave(
  filename = "output/figures/fig_sudan_panelB.pdf",
  plot = p_b,
  width = 8,
  height = 4,
  dpi = 300,
)

annote_state <- head(names(sort(colMeans(fit$effects$spill))), 3)

annote_spill <- create_spillover_annotations(
  fit = fit,
  time_point = 2012, # 注釈を打つ基準年
  units_to_label = annote_state,
  # Y座標が近い場合、ここで垂直方向に少しずらす
  vjust_list = c(0.2, 0.5, 0.8)
)

p_spill <- plot_spillover_panel(
  fit = fit,
  vline_at = 2011, # 垂直線
  annotations_df = annote_spill,
  skip_year = 2011,
  main_title = ""
)

p_spill

ggsave(
  filename = "output/figures/fig_sudan_panel_spillover.pdf",
  plot = p_spill,
  width = 8,
  height = 4,
  dpi = 300,
)


## ----fig2---------------------------------------------------------------------
# Robustness Check

# Prior Sensitivity
grid <- data.frame(
  a0 = c(3, 5, 2),
  b0 = c(1, 1, 0.5),
  rho_lo = c(-0.5, -0.3, -0.7),
  rho_hi = c(0.5, 0.3, 0.7),
  step_rho = c(0.05, 0.03, 0.07)
)

Yc_pre_obs <- panel_df %>%
  filter(
    (country != "Sudan") & (year < 2011)
  ) %>%
  pivot_wider(
    id_cols = year,
    names_from = country,
    values_from = colnames(panel_df)[4]
  ) %>%
  select(-year) %>%
  as.matrix()

Y0_pre <- panel_df %>%
  filter(
    (country == "Sudan") & (year < 2011)
  ) %>%
  select(colnames(panel_df)[4]) %>%
  as.matrix()

Xc_pre <- panel_df %>%
  filter(
    (country != "Sudan") & (year < 2011)
  ) %>%
  pivot_wider(
    id_cols = year,
    names_from = country,
    values_from = colnames(panel_df)[4]
  ) %>%
  select(-year) %>%
  as.matrix()

dim(Xc_pre) <- c(nrow(Xc_pre), ncol(Xc_pre), 1)

alpha_hat <- colMeans(fit$alpha_draws)

sens <- prior_sensitivity(
  Yc_obs = Yc_pre_obs,
  W_raw = W,
  w_raw = w,
  alpha_hat_scaled = alpha_hat,
  Xc_pre = Xc_pre,
  p = 0L,
  grid = grid,
  M_burn = sens_burn,
  M_keep = sens_keep
)

sens$theta_table

colnames(sens$theta_table) <- gsub("_", "\\\\_", colnames(sens$theta_table))

sens$theta_table %>%
  kable(
    "latex",
    booktabs = TRUE,
    caption = "Prior Sensitivity Analysis (Sudan Split)",
    label = "tab:sudan_sens",
    escape = FALSE
  ) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  save_kable("output/tables/sudan_sens_table.tex")

write.csv(
  sens$theta_table,
  "output/tables/sudan_sensitivity.csv",
  row.names = FALSE
)


## -----------------------------------------------------------------------------
# Prior Predictive

create_ppa_table <- function(ppc_results) {
  sim_stats <- ppc_results$stat

  obs_stats_vec <- as.numeric(ppc_results$observed)
  obs_names <- names(ppc_results$observed)

  if (is.null(obs_names)) {
    obs_names <- colnames(sim_stats)
  }

  p_values <- vapply(
    seq_along(obs_names),
    function(i) {
      stat_name <- obs_names[i]
      colnames(sim_stats) <- gsub("\\\\", "", colnames(sim_stats))
      sim_col <- sim_stats[[stat_name]]
      obs_val <- obs_stats_vec[i]

      sim_col_finite <- sim_col[is.finite(sim_col)]
      if (length(sim_col_finite) == 0 || !is.finite(obs_val)) {
        return(NA_real_)
      }

      mean(sim_col_finite <= obs_val)
    },
    numeric(1)
  )

  summary_table <- data.frame(
    Statistic = obs_names,
    Observed_Value_h_y_obs = obs_stats_vec,
    P_Value_P_h_y_obs = p_values,
    row.names = NULL
  )

  return(summary_table)
}

ppc <- prior_predictive(
  Y0_pre = Y0_pre,
  Yc_obs = Yc_pre_obs,
  W_raw = W,
  w_raw = w,
  alpha_hat_scaled = alpha_hat,
  Xc_pre = Xc_pre,
  p = 0L,
  a0 = 3,
  b0 = 1,
  rho_support = c(-0.99, 0.99),
  R = ppc_rep
)

stat_keys <- c(
  "yc_mean",
  "log_yc_var",
  "spatial_quadratic",
  "corr_y0_wyc",
  "ac1",
  "ac2",
  "pve_pc1",
  "avg_skewness",
  "avg_kurtosis"
)

resolve_stat_cols <- function(available_names, target_keys) {
  normalized <- gsub("\\\\", "", available_names)
  idx <- match(target_keys, normalized)
  if (any(is.na(idx))) {
    missing_keys <- target_keys[is.na(idx)]
    stop(sprintf(
      "Missing prior-predictive statistics: %s",
      paste(missing_keys, collapse = ", ")
    ))
  }
  available_names[idx]
}

sim_stat_cols <- resolve_stat_cols(colnames(ppc$stat), stat_keys)


sim_long <- ppc$stat %>%
  as.data.frame() %>%
  pivot_longer(
    cols = all_of(sim_stat_cols),
    names_to = "statistic",
    values_to = "simulated_value"
  ) %>%
  mutate(
    statistic = gsub("\\\\", "", statistic),
    statistic = gsub("_", " ", statistic)
  )

ppc_observed_tmp <- ppc$observed %>% t()
obs_stat_cols <- resolve_stat_cols(colnames(ppc_observed_tmp), stat_keys)

obs_long <- ppc_observed_tmp %>%
  as.data.frame() %>%
  pivot_longer(
    cols = all_of(obs_stat_cols),
    names_to = "statistic",
    values_to = "observed_value"
  ) %>%
  mutate(
    statistic = gsub("\\\\", "", statistic),
    statistic = gsub("_", " ", statistic)
  )

# --- 2. ggplot オブジェクトの作成 ---

# 'scales = "free_x"' が重要。これにより各プロットのX軸が独立します
gg <- ggplot(sim_long, aes(x = simulated_value)) +
  # 灰色のヒストグラム (事前予測分布)
  geom_histogram(
    bins = 40, # ビンの数を調整
    fill = "grey90",
    color = "white",
    alpha = 0.8
  ) +

  # 赤色の縦線 (観測値)
  # obs_long を geom_vline の data として指定
  geom_vline(
    data = obs_long,
    aes(xintercept = observed_value),
    color = "red",
    linewidth = 1.0 # 線の太さを調整
  ) +

  # 統計量ごとに 3x3 のグリッドでファセット化
  facet_wrap(
    ~statistic,
    scales = "free", # X軸もY軸も統計量ごとにスケールを自動調整
    ncol = 3 # 3列で表示
  ) +

  # 論文用のテーマ
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold") # 各パネルのタイトル (e.g., "yc_mean")
  ) +
  labs(
    x = "Statistic Value",
    y = "Frequency"
  )

print(gg)

ggsave(
  filename = "output/figures/fig_sudan_prior_predictive.pdf",
  plot = gg,
  width = 8.0,
  height = 4.0,
  dpi = 300,
)

ppc_table <- create_ppa_table(ppc)

rename_vec <- c(
  yc_mean = "Mean",
  log_yc_var = "log-Var",
  spatial_quadratic = "Spatial Quad.",
  corr_y0_wyc = "Corr($Y_0$, $w'Y$)",
  ac1 = "AC(1)",
  ac2 = "AC(2)",
  pve_pc1 = "PVE(PC1)",
  avg_skewness = "Skewness",
  avg_kurtosis = "Kurtosis"
)
ppc_table$Statistic <- rename_vec[ppc_table$Statistic]
ppc_table$Observed_Value_h_y_obs <- round(ppc_table$Observed_Value_h_y_obs, 3)
ppc_table$P_Value_P_h_y_obs <- round(ppc_table$P_Value_P_h_y_obs, 3)

print(ppc_table)

write.csv(
  ppc_table,
  "output/tables/sudan_ppa_summary_table.csv",
  row.names = FALSE
)

latex_col_names <- c("Statistic", "Observed", "$P(h \\leq h(y^{o}) | A)$")

ppc_table <- ppc_table %>% mutate(Statistic = gsub("_", "\\\\_", Statistic))

ppa_kable_object <- ppc_table %>%
  # Round all numeric columns to 3 decimal places for a clean look
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%

  kable(
    "latex",
    booktabs = TRUE,
    caption = "Prior Predictive Analysis Summary Table (Sudan Split)", # Updated caption
    label = "tab:sudan_ppa", # Updated label
    escape = FALSE, # IMPORTANT: Renders LaTeX in col.names
    col.names = latex_col_names, # Use LaTeX column names
    align = "lrr" # Align columns (Left, Right, Right)
  ) %>%
  kable_styling(latex_options = c("hold_position"))

# Save the table to a .tex file
save_kable(ppa_kable_object, "output/tables/sudan_ppa_summary_table.tex")

# (Optional) Print the LaTeX code to the console to check
print(ppa_kable_object)


## -----------------------------------------------------------------------------
# ---- packages ----
library(coda)
library(mcmcse)
library(posterior)

# ---- helpers ----
.check_draws_matrix <- function(x) {
  # Ensure a numeric matrix with column names
  stopifnot(is.matrix(x), is.numeric(x))
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("param[", seq_len(ncol(x)), "]")
  }
  x
}

# ============================================
# Single-chain diagnostics
# Input: matrix (iterations x parameters)
# Output: data.frame of summaries; attributes contain classical tests
# ============================================
mcmc_diagnostics_single <- function(draws_matrix, label = "chain1") {
  # draws_matrix: numeric matrix (iterations x parameters)
  Z <- .check_draws_matrix(draws_matrix)
  mmc <- coda::mcmc(Z)

  # means, sd
  mu <- colMeans(Z)
  sds <- apply(Z, 2, sd)

  # ESS (coda) and lag-1 ACF
  ess_coda <- coda::effectiveSize(mmc)
  acf_lag1 <- apply(Z, 2, function(z) {
    acf(z, plot = FALSE, lag.max = 1)$acf[2L]
  })

  # MCSE (mcmcse, batch means)
  mcse_mean <- sapply(seq_len(ncol(Z)), function(j) mcmcse::mcse(Z[, j])$se)
  ess_bm <- sapply(seq_len(ncol(Z)), function(j) mcmcse::ess(Z[, j]))

  # posterior::ess_bulk / ess_tail; Rhat is NA for single chain
  ds <- posterior::as_draws_matrix(Z)
  ess_bulk <- posterior::ess_bulk(ds)
  ess_tail <- posterior::ess_tail(ds)
  rhat <- posterior::rhat(ds) # single chain -> NA

  # Classical diagnostics (stored in attributes)
  gz <- try(suppressWarnings(coda::geweke.diag(mmc)$z), silent = TRUE)
  hd <- try(suppressWarnings(coda::heidel.diag(mmc)), silent = TRUE)
  rft <- try(suppressWarnings(coda::raftery.diag(mmc)), silent = TRUE)

  out <- data.frame(
    param = colnames(Z),
    mean = as.numeric(mu),
    sd = as.numeric(sds),
    ess_coda = as.numeric(ess_coda),
    ess_bulk = as.numeric(ess_bulk),
    ess_tail = as.numeric(ess_tail),
    ess_bm = as.numeric(ess_bm),
    mcse_mean = as.numeric(mcse_mean),
    acf_lag1 = as.numeric(acf_lag1),
    rhat = as.numeric(rhat)
  )
  attr(out, "gewekeZ") <- gz
  attr(out, "heidel") <- hd
  attr(out, "raftery") <- rft
  attr(out, "label") <- label
  out
}


# ============================================
# Optional: quick plots if bayesplot is available
# Input: draws (matrix) or list of matrices
# ============================================
mcmc_quick_plots <- function(draws, params = NULL, max_params = 8) {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    message("bayesplot not installed; skipping plots.")
    return(invisible(NULL))
  }
  pick_cols <- function(m) {
    idx <- seq_len(ncol(m))
    if (!is.null(params)) {
      idx <- which(colnames(m) %in% params)
    }
    idx <- head(idx, max_params)
    m[, idx, drop = FALSE]
  }

  if (is.matrix(draws)) {
    Z <- pick_cols(.check_draws_matrix(draws))
    bayesplot::mcmc_trace(Z)
    bayesplot::mcmc_acf(Z)
  } else if (is.list(draws)) {
    L <- lapply(draws, function(m) pick_cols(.check_draws_matrix(m)))
    # to array: iterations x parameters x chains
    iters <- min(sapply(L, nrow))
    P <- ncol(L[[1L]])
    arr <- array(NA_real_, dim = c(iters, P, length(L)))
    for (c in seq_along(L)) {
      arr[,, c] <- L[[c]][seq_len(iters), , drop = FALSE]
    }
    colnames(arr) <- colnames(L[[1L]])
    dm <- posterior::as_draws_array(arr)
    bayesplot::mcmc_trace(dm)
    bayesplot::mcmc_acf(dm)
    bayesplot::mcmc_rhat(dm)
  } else {
    stop("draws must be a matrix or a list of matrices.")
  }
  invisible(NULL)
}

D <- cbind(fit$alpha_draws, fit$rho_draws)
diag1 <- mcmc_diagnostics_single(D, label = "posterior_run")
diag1


## -----------------------------------------------------------------------------
w <- as.numeric(fit$inputs$w)
country_names <- fit$inputs$units$controls

w_df <- data.frame(
  country = country_names,
  weight = w
)

w_df_filtered <- w_df %>%
  # filter(weight > 1e-8) %>%
  arrange(weight)

gg_w_bar <- ggplot(
  w_df_filtered,
  aes(
    x = weight,
    y = reorder(country, weight)
  )
) +
  geom_col(fill = "gray40", width = 0.7) +
  labs(
    x = "Weight",
    y = "Control country"
  ) +
  theme_minimal() +
  xlim(0.0, NA)

print(gg_w_bar)

ggsave(
  "output/figures/fig_sudan_trade.pdf",
  height = 8,
  width = 6
)


## -----------------------------------------------------------------------------

proposed_w <- fit$alpha_hat
countries <- names(proposed_w)

scm_vec <- compute_scm_weights(Y0_pre, Yc_pre_obs)
names(scm_vec) <- countries

df_weights <- tibble(
  country = countries,
  Proposed = as.numeric(proposed_w[countries]),
  SCM = as.numeric(scm_vec[countries])
) %>%
  pivot_longer(
    cols = c(Proposed, SCM),
    names_to = "method",
    values_to = "weight"
  ) %>%
  mutate(
    country = factor(country, levels = sort(unique(country)))
  )

ggplot(df_weights, aes(x = weight, y = country, fill = method)) +
  geom_col_pattern(
    aes(pattern = method), # パターンだけを method に割り当て
    position = position_dodge(width = 0.7),
    colour = NA,
    pattern_fill = "black",
    pattern_colour = NA,
    pattern_angle = 30,
    pattern_density = 0.1,
    pattern_spacing = 0.01
  ) +
  scale_pattern_manual(
    values = c("SCM" = "none", "Proposed" = "stripe")
  ) +
  coord_flip() +
  labs(x = NULL, y = "Country", fill = "method") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  "output/figures/fig_sudan_weight.pdf",
  height = 5.5,
  width = 11
)
