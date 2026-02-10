## ----setup, message=FALSE, warning=FALSE--------------------------------------
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(knitr)
library(kableExtra)

run_mode <- tolower(Sys.getenv("SCSPILL_MODE", "full"))
if (!run_mode %in% c("full", "smoke")) {
  stop("SCSPILL_MODE must be either 'full' or 'smoke'.")
}
if (run_mode == "smoke") {
  Ns <- c(16)
  T0_vals <- c(20)
  T1_vals <- c(20)
  rho_vals <- c(-0.3, 0.0, 0.3)
  sims_per <- 20L
  mcmc_iter <- 800L
  burn_iter <- 200L
} else {
  Ns <- c(16, 36, 64)
  T0_vals <- c(20, 50)
  T1_vals <- c(20, 50)
  rho_vals <- c(-0.8, -0.3, -0.1, 0.0, 0.1, 0.3, 0.8)
  sims_per <- 1000L
  mcmc_iter <- 2500L
  burn_iter <- 500L
}

dir.create("output/tables/mc_result", recursive = TRUE, showWarnings = FALSE)

source("R/01_utils.R")
source("R/10_sc_spillover.R")
source("R/21_mcmc_alpha.R")
source("R/22_mcmc_sar.R")
source("R/30_simulation_.R")
source("R/40_geweke_latest.R")
source("R/41_robustness_check.R")
Rcpp::sourceCpp("src/20_mcmc.cpp")
Rcpp::sourceCpp("src/40_geweke_latest.cpp")

# Assume the simulation engine functions are already part of the package
set.seed(20251030)
message(sprintf("[seed-init] script=simulation_main global_seed=%d mode=%s", 20251030L, run_mode))


## -----------------------------------------------------------------------------
make_rook_W <- function(nrow, ncol, normalize = FALSE) {
  stopifnot(nrow >= 1L, ncol >= 1L)
  idx <- matrix(seq_len(nrow * ncol), nrow, ncol, byrow = TRUE)
  N <- nrow * ncol
  W <- matrix(0, N, N)

  for (r in 1:nrow) {
    for (c in 1:ncol) {
      i <- idx[r, c]
      if (r > 1) {
        W[i, idx[r - 1, c]] <- 1
      }
      if (r < nrow) {
        W[i, idx[r + 1, c]] <- 1
      }
      if (c > 1) {
        W[i, idx[r, c - 1]] <- 1
      }
      if (c < ncol) W[i, idx[r, c + 1]] <- 1
    }
  }

  if (normalize) {
    rs <- rowSums(W)
    W <- W / pmax(rs, 1) # 0割回避
  }
  W
}

run_scenario_many <- function(
  n_sims,
  dgp_args,
  seeds = NULL,
  M = 5000,
  burn = 2000,
  step_rho = 0.02
) {
  res <- run_many_sim(
    n_sims = n_sims,
    dgp_args = dgp_args,
    seeds = seeds,
    M = M,
    burn = burn,
    step_rho = step_rho
  )
  list(
    summary = summarize_many(res),
    raw = res
  )
}

mc_grid_study <- function(
  Ns,
  T0s,
  rhos,
  T1,
  sims_per = 100,
  K = 1,
  beta = c(1.0),
  sigma2 = 1.0,
  treated_idx = 1:4,
  M = 5000,
  burn = 2000,
  step_rho = 0.02,
  seeds = NULL
) {
  # グリッド展開
  grid <- expand.grid(N = Ns, T0 = T0s, rho = rhos, stringsAsFactors = FALSE)

  # 乱数種
  if (is.null(seeds)) {
    seeds <- sample.int(.Machine$integer.max, nrow(grid) * sims_per)
  }
  # シナリオごとの種をスライス
  get_seeds_i <- function(i) {
    idx <- ((i - 1L) * sims_per + 1L):(i * sims_per)
    seeds[idx]
  }

  # 実行
  details <- vector("list", nrow(grid))
  summaries <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    N <- grid$N[i]
    T0 <- grid$T0[i]
    rho <- grid$rho[i]

    # N は正方格子前提（論文の体裁を踏襲）。
    m <- round(sqrt(N))
    if (m * m != N) {
      stop(sprintf(
        "N=%d は完全平方数ではありません（%dx%dにできません）。",
        N,
        m,
        m
      ))
    }

    # W と w
    W <- make_rook_W(m, m)
    w <- numeric(N)
    w[intersect(treated_idx, seq_len(N))] <- 1

    # 真の α
    alpha_true <- rep(0, N)
    alpha_true[1] <- 0.5
    alpha_true[2] <- -0.2
    alpha_true[3:4] <- 0.4
    alpha_true[5:10] <- 0.1 / 6

    # DGP 引数
    dgp_args <- list(
      T0 = T0,
      T1 = T1,
      N = N,
      W = W,
      w = w,
      rho = rho,
      sigma2 = sigma2,
      alpha = alpha_true,
      K = K,
      beta = beta
    )

    # 多回実行
    process_time <- system.time({
      out <- run_scenario_many(
        n_sims = sims_per,
        dgp_args = dgp_args,
        seeds = get_seeds_i(i),
        M = M,
        burn = burn,
        step_rho = step_rho
      )
    })
    print(process_time)

    # シナリオ情報を列として付与
    smry <- out$summary
    smry$N <- N
    smry$T0 <- T0
    smry$T1 <- T1
    smry$rho <- rho
    # 列順を調整（読みやすさ）
    smry <- smry[, c(
      "N",
      "T0",
      "T1",
      "rho",
      "method",
      "bias_point",
      "rmse_point",
      "cover95_point"
    )]

    details[[i]] <- out$raw
    summaries[[i]] <- smry
  }

  list(
    summary = do.call(rbind, summaries),
    details = details,
    grid = grid
  )
}


## ----run-grid-----------------------------------------------------------------
output_dir <- "output/tables/mc_result"

for (N in Ns) {
  for (T0 in T0_vals) {
    for (T1 in T1_vals) {
      study <- mc_grid_study(
        Ns = c(N), # grid size: 4×4 and 6×6
        T0s = c(T0), # pre-period length
        rhos = rho_vals,
        T1 = T1, # post-treatment length (20 or 50)
        sims_per = sims_per,
        K = 1, # number of covariates
        beta = c(1.0), # true regression slope
        sigma2 = 0.1, # innovation variance
        treated_idx = 1:4, # treated units
        M = mcmc_iter,
        burn = burn_iter,
        step_rho = 0.05 # MH step size
      )
      output_path <- paste0(
        output_dir,
        "/mc_study_N=",
        N,
        "T0=",
        T0,
        "T1=",
        T1,
        ".csv"
      )
      write.csv(study$summary, output_path, row.names = FALSE)
      rm(study)
    }
  }
}



## -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(glue)


read_mc_results <- function(dir = "output/tables/mc_result") {
  files <- list.files(
    dir,
    pattern = "^mc_study_N=\\d+T0=\\d+T1=\\d+\\.csv$",
    full.names = TRUE
  )

  meta <- tibble(file = files) %>%
    mutate(
      fname = basename(file),
      N = as.integer(str_match(fname, "N=(\\d+)")[, 2]),
      T0 = as.integer(str_match(fname, "T0=(\\d+)")[, 2]),
      T1 = as.integer(str_match(fname, "T1=(\\d+)")[, 2]),
      data = map(
        file,
        ~ readr::read_csv(.x, show_col_types = FALSE) %>%
          select(-any_of(c("N", "T0", "T1")))
      )
    )
  res <- meta %>%
    unnest(cols = data) %>%
    rename(
      Bias = bias_point,
      RMSE = rmse_point,
      CoverageRate = cover95_point
    )
  res <- res %>%
    mutate(method = ifelse(method == "SCSPILL", "Proposed", method))
  res
}

mc_res <- read_mc_results("output/tables/mc_result")

make_table1_panel_lines <- function(res, T0_val, T1_val, panel_label) {
  df <- res %>%
    filter(T0 == T0_val, T1 == T1_val) %>%
    filter(method %in% c("Proposed", "SCM", "BSCM")) %>%
    mutate(
      method = factor(method, levels = c("Proposed", "SCM", "BSCM")),
      rho = as.numeric(rho)
    )

  N_vals <- sort(unique(df$N))
  rho_vals <- sort(unique(df$rho))
  rho_vals <- round(rho_vals, digits = 1)

  fmt <- function(x) formatC(x, format = "f", digits = 3)

  col_spec <- paste0("ll", paste(rep("rrr", length(N_vals)), collapse = ""))

  cmid_lines <- map_chr(seq_along(N_vals), function(i) {
    start <- 3 + 3 * (i - 1)
    end <- start + 2
    glue("\\cmidrule(lr){{{start}-{end}}}")
  }) %>%
    paste(collapse = "")

  lines <- c(
    glue("\\caption*{{Panel ({panel_label}): $T_0={T0_val}$}}"),
    glue("\\label{{5.2 Simulation panel {panel_label}}}"),
    glue("\\begin{{tabular}}{{{col_spec}}}"),
    "\\hline",
    "\\multicolumn{2}{l}{} & ",
    paste(
      map_chr(N_vals, ~ glue("\\multicolumn{{3}}{{c}}{{$N={.x}$}}")),
      collapse = " & "
    ),
    " \\\\",
    cmid_lines,
    "\n",
    "\\multicolumn{1}{l}{} & \\multicolumn{1}{l}{$\\rho$} & ",
    paste(rep(c("Proposed", "SCM", "BSCM"), length(N_vals)), collapse = " & "),
    " \\\\ \\hline"
  )

  first <- TRUE
  for (r in rho_vals) {
    if (first) {
      line <- " & "
      first <- FALSE
    } else {
      if (r == 0.0) {
        line <- "Bias & "
      } else {
        line <- "& "
      }
    }
    line <- paste0(line, round(r, digits = 1))
    for (N0 in N_vals) {
      tmp <- df %>%
        filter(N == N0, rho == r) %>%
        arrange(method)
      vals <- fmt(tmp$Bias)
      line <- paste0(line, " & ", paste(vals, collapse = " & "))
    }
    lines <- c(lines, paste0(line, " \\\\"))
  }

  lines <- c(lines, "\\hdashline")

  first <- TRUE
  for (r in rho_vals) {
    if (first) {
      line <- " & "
      first <- FALSE
    } else {
      if (r == 0.0) {
        line <- "RMSE & "
      } else {
        line <- "& "
      }
    }
    line <- paste0(line, round(r, digits = 1))
    for (N0 in N_vals) {
      tmp <- df %>%
        filter(N == N0, rho == r) %>%
        arrange(method)
      vals <- fmt(tmp$RMSE)
      line <- paste0(line, " & ", paste(vals, collapse = " & "))
    }
    lines <- c(lines, paste0(line, " \\\\"))
  }

  lines <- c(lines, "\\hline", "\\end{tabular}")

  lines
}

make_table2_lines <- function(mc_res, method_proposed = "Proposed") {
  mc_res <- mc_res %>%
    mutate(rho = as.numeric(rho))

  cov_df <- mc_res %>%
    filter(method == method_proposed) %>%
    group_by(T0, N, rho) %>%
    summarise(coverage = mean(CoverageRate, na.rm = TRUE), .groups = "drop")
  t0_vals <- sort(unique(cov_df$T0))
  n_vals <- sort(unique(cov_df$N))
  rho_vals <- sort(unique(cov_df$rho))
  block_width <- length(n_vals)

  col_spec <- paste0("l", paste(rep("r", block_width * length(t0_vals)), collapse = ""))
  header_top <- paste(
    purrr::map_chr(t0_vals, ~ sprintf("\\multicolumn{%d}{c}{$T_{0}=%d$}", block_width, .x)),
    collapse = " & "
  )
  cmid <- paste(
    purrr::map_chr(seq_along(t0_vals), function(i) {
      start <- 2 + (i - 1) * block_width
      end <- start + block_width - 1
      sprintf("\\cmidrule(lr){%d-%d}", start, end)
    }),
    collapse = " "
  )
  header_bottom <- paste(rep(sprintf("$N=%d$", n_vals), length(t0_vals)), collapse = " & ")

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Coverage Rate of 95\\% Credible Interval for the Proposed Method}",
    "\\label{5.2 Simulation coverage rate}",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\hline",
    paste("$\\rho$ &", header_top, "\\\\"),
    cmid,
    paste("$\\rho$ &", header_bottom, "\\\\ \\hline")
  )

  for (r in rho_vals) {
    vals <- character(0)
    for (t0 in t0_vals) {
      for (n0 in n_vals) {
        v <- cov_df %>%
          filter(T0 == t0, N == n0, rho == r) %>%
          pull(coverage)
        vals <- c(vals, if (length(v) == 0) "" else sprintf("%.3f", v[1]))
      }
    }
    lines <- c(lines, paste(sprintf("% .1f", r), "&", paste(vals, collapse = " & "), "\\\\"))
  }

  notes <- c(
    "\\hline",
    "\\end{tabular}",
    "\\begin{tablenotes}",
    "{\\footnotesize",
    "\\item Notes: This table shows the coverage rate of $95$\\% credible interval for the proposed method for each scenario of $\\rho$, $T_0$, and  $N$. The coverage rate is computed over 1000 simulations.",
    "}",
    "\\end{tablenotes}",
    "\\end{table}"
  )

  c(lines, notes)
}

write_simulation_tables <- function(
  mc_res,
  T1_val = 20,
  file = "output/tables/simulation_tables.tex"
) {
  available_t0 <- sort(unique(mc_res$T0))
  lines_panel_a <- make_table1_panel_lines(
    res = mc_res,
    T0_val = available_t0[1],
    T1_val = T1_val,
    panel_label = "a"
  )

  lines_panel_b <- character(0)
  if (length(available_t0) >= 2) {
    lines_panel_b <- c(
      "",
      "\\bigskip",
      "",
      make_table1_panel_lines(
        res = mc_res,
        T0_val = available_t0[2],
        T1_val = T1_val,
        panel_label = "b"
      )
    )
  }

  # Table 1 全体
  table1_lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Simulation Results for Bias and RMSE}",
    "\\label{5.2 Simulation}",
    lines_panel_a,
    lines_panel_b,
    "",
    "\\medskip",
    "\\begin{tablenotes}",
    "{\\footnotesize",
    "\\item Notes: Panels (a) and (b) show the simulation results for the bias and RMSE for $T_0=20$ and $50$, respectively. Each panel shows the simulation results for each of the proposed method, SCM, and BSCM, and each of $\\rho \\in \\{-0.8,-0.3, -0.1, 0.0, 0.1, 0.3,0.8\\}$ and $N \\in \\{16,36,64\\}$.",
    "}",
    "\\end{tablenotes}",
    "\\end{table}"
  )

  table2_lines <- make_table2_lines(mc_res, method_proposed = "Proposed")

  all_lines <- c(table1_lines, "", table2_lines)
  writeLines(all_lines, con = file)

  invisible(file)
}

write_simulation_tables(
  mc_res,
  T1_val = 20,
  file = "output/tables/simulation_results.tex"
)
