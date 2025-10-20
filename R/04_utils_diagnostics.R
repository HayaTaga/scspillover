#' scspill: Diagnostics (trace plots, multi-parameter)
#'
#' @param object scspill オブジェクト
#' @param what   "trace" のみをサポート（将来拡張用）
#' @param which_alpha  alpha の表示ユニット（文字=ユニット名 or 数値=列番号）。NULL なら top_n_alpha を自動選抜
#' @param top_n_alpha  which_alpha 未指定時に |alpha_hat| の大きい順に抽出する数（既定 6）
#' @param which_beta   beta の表示（colnames or 列番号）。NULL なら top_n_beta を自動選抜（あれば）
#' @param top_n_beta   which_beta 未指定時の抽出数（既定 6）
#' @export
diagnostics.scspill <- function(
  object,
  what = c("trace"),
  which_alpha = NULL,
  top_n_alpha = 6,
  which_beta = NULL,
  top_n_beta = 6
) {
  what <- match.arg(what)
  if (what != "trace") {
    stop("Currently only what='trace' is supported.")
  }

  # ---- 収集：rho / alpha / sigma2 / tau2 / beta（存在すれば）を縦持ちで結合 ----
  out_list <- list()

  # rho
  if (!is.null(object$rho_draws)) {
    out_list[["rho"]] <- data.frame(
      iter = seq_along(object$rho_draws),
      value = as.numeric(object$rho_draws),
      series = "rho"
    )
  }

  # alpha（上位 or 指定）
  if (!is.null(object$alpha_draws)) {
    unit_names <- {
      uc <- tryCatch(object$inputs$units$control, error = function(e) NULL)
      if (!is.null(uc)) {
        as.character(uc)
      } else {
        paste0("unit_", seq_len(ncol(object$alpha_draws)))
      }
    }
    # 選抜
    if (is.null(which_alpha)) {
      ord <- order(-abs(as.numeric(object$alpha_hat)))
      sel <- head(ord, n = min(top_n_alpha, length(ord)))
    } else if (is.character(which_alpha)) {
      m <- match(which_alpha, unit_names)
      if (anyNA(m)) {
        stop(
          "Unknown unit(s) in 'which_alpha': ",
          paste(which_alpha[is.na(m)], collapse = ", ")
        )
      }
      sel <- as.integer(m)
    } else if (is.numeric(which_alpha)) {
      sel <- as.integer(which_alpha)
      if (any(sel < 1 | sel > ncol(object$alpha_draws))) {
        stop("'which_alpha' indices out of range.")
      }
    } else {
      stop("'which_alpha' must be NULL, character, or numeric.")
    }
    # ロング化
    iters <- seq_len(nrow(object$alpha_draws))
    for (j in sel) {
      out_list[[paste0("alpha[", unit_names[j], "]")]] <- data.frame(
        iter = iters,
        value = as.numeric(object$alpha_draws[, j]),
        series = paste0("alpha[", unit_names[j], "]")
      )
    }
  }

  # sigma2 / tau2（あれば sar に格納していることを想定）
  if (!is.null(object$sar) && !is.null(object$sar$sigma2_draws)) {
    out_list[["sigma2"]] <- data.frame(
      iter = seq_along(object$sar$sigma2_draws),
      value = as.numeric(object$sar$sigma2_draws),
      series = "sigma2"
    )
  }
  if (!is.null(object$sar) && !is.null(object$sar$tau2_draws)) {
    out_list[["tau2"]] <- data.frame(
      iter = seq_along(object$sar$tau2_draws),
      value = as.numeric(object$sar$tau2_draws),
      series = "tau2"
    )
  }

  # beta（あれば; M x K を想定）
  if (!is.null(object$sar) && !is.null(object$sar$beta)) {
    beta_draws <- object$sar$beta
    K <- ncol(beta_draws)
    beta_names <- colnames(beta_draws)
    if (is.null(beta_names)) {
      beta_names <- paste0("beta_", seq_len(K))
    }

    if (is.null(which_beta)) {
      # 事後平均の |.| 大きい順に抽出
      bm <- colMeans(beta_draws)
      selb <- order(-abs(bm))
      selb <- head(selb, n = min(top_n_beta, length(selb)))
    } else if (is.character(which_beta)) {
      mb <- match(which_beta, beta_names)
      if (anyNA(mb)) {
        stop(
          "Unknown name(s) in 'which_beta': ",
          paste(which_beta[is.na(mb)], collapse = ", ")
        )
      }
      selb <- as.integer(mb)
    } else if (is.numeric(which_beta)) {
      selb <- as.integer(which_beta)
      if (any(selb < 1 | selb > K)) stop("'which_beta' indices out of range.")
    } else {
      stop("'which_beta' must be NULL, character, or numeric.")
    }

    iters <- seq_len(nrow(beta_draws))
    for (k in selb) {
      out_list[[paste0("beta[", beta_names[k], "]")]] <- data.frame(
        iter = iters,
        value = as.numeric(beta_draws[, k]),
        series = paste0("beta[", beta_names[k], "]")
      )
    }
  }

  if (length(out_list) == 0L) {
    stop("No traceable parameters found in object.")
  }

  df <- do.call(rbind, out_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(iter, value)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~series, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Iteration",
      y = "Value",
      title = "Trace plots (rho / alpha / sigma2 / tau2 / beta)"
    )
  return(p)
}
