# ---------------- 共通ユーティリティ ----------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

.get_times_pre <- function(fit) {
  # 優先: inputs$times_pre, 次点: rownames(Yc_pre), 最後: 1:T0
  if (!is.null(fit$inputs$times_pre)) {
    return(as.vector(fit$inputs$times_pre))
  }
  if (!is.null(rownames(fit$inputs$Yc_pre))) {
    rn <- rownames(fit$inputs$Yc_pre)
    if (suppressWarnings(all(!is.na(as.numeric(rn))))) {
      return(as.numeric(rn))
    }
    return(rn)
  }
  seq_len(nrow(fit$inputs$Yc_pre))
}

.get_times_post <- function(fit) {
  if (!is.null(fit$inputs$times_post)) {
    return(as.vector(fit$inputs$times_post))
  }
  if (!is.null(rownames(fit$inputs$Yc_post))) {
    rn <- rownames(fit$inputs$Yc_post)
    if (suppressWarnings(all(!is.na(as.numeric(rn))))) {
      return(as.numeric(rn))
    }
    return(rn)
  }
  seq_len(nrow(fit$inputs$Yc_post))
}

.get_units_control <- function(fit, n_cols) {
  if (!is.null(fit$inputs$units) && is.list(fit$inputs$units)) {
    uc <- fit$inputs$units$control
    if (!is.null(uc)) return(as.character(uc))
  }
  cn <- colnames(fit$inputs$Yc_post) %||% colnames(fit$inputs$Yc_pre)
  if (!is.null(cn)) {
    return(as.character(cn))
  }
  paste0("unit_", seq_len(n_cols))
}

# -------------- counterfactual（事後平均/区間を計算） ----------------
# ※ posterior_effects() の同じ同定式を用い、treated の反事実 y0_cf を再構成
#   Ainv = (I_N - r (w a' + W))^{-1}
#   tmp  = Ainv { (I_N - r W) y_c - r w y0 }
#   effect = y0 - a' tmp  =>  counterfactual = a' tmp
#
# 返り値: data.frame(time, period, y_obs, y_cf_mean, y_cf_lo, y_cf_hi)

#' @keywords internal
scspill_counterfactual <- function(fit, cred = 0.95) {
  stopifnot(inherits(fit, "scspill"))

  Y0_pre <- as.numeric(fit$inputs$Y0_pre)
  Yc_pre <- as.matrix(fit$inputs$Yc_pre)
  storage.mode(Yc_pre) <- "double"
  Y0_post <- as.numeric(fit$inputs$Y0_post)
  Yc_post <- as.matrix(fit$inputs$Yc_post)
  storage.mode(Yc_post) <- "double"

  N <- ncol(Yc_pre)
  T0 <- nrow(Yc_pre)
  T1 <- nrow(Yc_post)

  w <- as.matrix(fit$inputs$w)
  storage.mode(w) <- "double"
  W <- as.matrix(fit$inputs$W)
  storage.mode(W) <- "double"

  alpha_draws <- as.matrix(fit$alpha_draws) # M x N
  rho_draws <- as.numeric(fit$rho_draws)
  M <- nrow(alpha_draws)

  IN <- diag(N)

  times_pre <- fit$inputs$times_pre
  times_post <- fit$inputs$times_post

  # ---- pre 期間の counterfactual（事後サンプル）----
  ycf_pre_draws <- matrix(NA_real_, nrow = T0, ncol = M)
  for (m in seq_len(M)) {
    a <- alpha_draws[m, ]
    r <- rho_draws[m]
    Ainv <- solve(IN - r * (w %*% t(a) + W))
    B <- (IN - r * W)
    for (t in seq_len(T0)) {
      yc <- Yc_pre[t, ]
      y0 <- Y0_pre[t]
      tmp <- Ainv %*% (B %*% yc - r * w * y0)
      ycf_pre_draws[t, m] <- as.numeric(crossprod(a, tmp))
    }
  }

  # ---- post 期間の counterfactual（事後サンプル）----
  ycf_post_draws <- matrix(NA_real_, nrow = T1, ncol = M)
  for (m in seq_len(M)) {
    a <- alpha_draws[m, ]
    r <- rho_draws[m]
    Ainv <- solve(IN - r * (w %*% t(a) + W))
    B <- (IN - r * W)
    for (t in seq_len(T1)) {
      yc <- Yc_post[t, ]
      y0 <- Y0_post[t]
      tmp <- Ainv %*% (B %*% yc - r * w * y0)
      ycf_post_draws[t, m] <- as.numeric(crossprod(a, tmp))
    }
  }

  lo_q <- (1 - cred) / 2
  hi_q <- 1 - lo_q
  # pre
  ycf_pre_mean <- rowMeans(ycf_pre_draws)
  ycf_pre_lo <- apply(ycf_pre_draws, 1, stats::quantile, probs = lo_q)
  ycf_pre_hi <- apply(ycf_pre_draws, 1, stats::quantile, probs = hi_q)
  # post
  ycf_post_mean <- rowMeans(ycf_post_draws)
  ycf_post_lo <- apply(ycf_post_draws, 1, stats::quantile, probs = lo_q)
  ycf_post_hi <- apply(ycf_post_draws, 1, stats::quantile, probs = hi_q)

  df_pre <- data.frame(
    time = times_pre,
    t_idx = seq_len(T0),
    period = "pre",
    y_obs = Y0_pre,
    y_cf_mean = ycf_pre_mean,
    y_cf_lo = ycf_pre_lo,
    y_cf_hi = ycf_pre_hi
  )
  df_post <- data.frame(
    time = times_post,
    t_idx = T0 + seq_len(T1),
    period = "post",
    y_obs = Y0_post,
    y_cf_mean = ycf_post_mean,
    y_cf_lo = ycf_post_lo,
    y_cf_hi = ycf_post_hi
  )
  rbind(df_pre, df_post)
}

# ---------------- tidy 化（従来：effects も使う） ----------------

#' @keywords internal
tidy_scspill <- function(fit) {
  stopifnot(inherits(fit, "scspill"))

  times_post <- .get_times_post(fit)
  units_ctrl <- .get_units_control(fit, ncol(fit$effects$spill$mean))

  # treated effect（post）
  df_treat <- data.frame(
    time = times_post,
    mean = fit$effects$treat$mean,
    lo = fit$effects$treat$lo,
    hi = fit$effects$treat$hi
  )

  # spill（post）
  sm <- fit$effects$spill$mean
  slo <- fit$effects$spill$lo
  shi <- fit$effects$spill$hi
  df_spill <- do.call(
    rbind,
    lapply(seq_len(ncol(sm)), function(j) {
      data.frame(
        time = times_post,
        unit = units_ctrl[j],
        mean = sm[, j],
        lo = slo[, j],
        hi = shi[, j]
      )
    })
  )

  # パラメータ要約
  df_param <- data.frame(
    param = "rho",
    mean = mean(fit$rho_draws),
    sd = stats::sd(fit$rho_draws),
    q025 = stats::quantile(fit$rho_draws, 0.025),
    q975 = stats::quantile(fit$rho_draws, 0.975)
  )

  weights <- data.frame(
    unit = units_ctrl,
    alpha = as.numeric(fit$alpha_hat)
  )

  list(treat = df_treat, spill = df_spill, params = df_param, weights = weights)
}

# ---------------- S3: plot.scspill（既定=前後含む“実測vs反事実”） ----------------

#' scspill の可視化
#'
#' @param x scspill オブジェクト
#' @param type "full"（既定: 実測vs反事実, 前後一括）, "effect"（postの効果推移）,
#'             "spill_top", "weights", "rho", "beta", "trace"
#' @param top_n spill_top の表示ユニット数
#' @param cred  信用水準（反事実リボン）
#' @param ...   予備
#' @export
plot.scspill <- function(
  x,
  type = c("full", "effect", "spill_top", "weights", "rho", "beta", "trace"),
  top_n = 8,
  cred = 0.95,
  ...
) {
  type <- match.arg(type)
  gg <- NULL

  if (type == "full") {
    cf <- scspill_counterfactual(x, cred = cred)

    df_pre <- subset(cf, period == "pre")
    df_post <- subset(cf, period == "post")

    gg <- ggplot2::ggplot(cf, ggplot2::aes(x = t_idx)) +
      # 観測（pre+post）— 実線
      ggplot2::geom_line(ggplot2::aes(y = y_obs), linewidth = 0.5) +
      # 反事実（pre）— リボン + 破線
      ggplot2::geom_ribbon(
        data = df_pre,
        ggplot2::aes(ymin = y_cf_lo, ymax = y_cf_hi),
        alpha = 0.12
      ) +
      ggplot2::geom_line(
        data = df_pre,
        ggplot2::aes(y = y_cf_mean),
        linetype = "22"
      ) +
      # 反事実（post）— リボン + 破線
      ggplot2::geom_ribbon(
        data = df_post,
        ggplot2::aes(ymin = y_cf_lo, ymax = y_cf_hi),
        alpha = 0.18
      ) +
      ggplot2::geom_line(
        data = df_post,
        ggplot2::aes(y = y_cf_mean),
        linetype = "22"
      ) +
      # 介入境界（t_idx = T0）
      ggplot2::geom_vline(xintercept = max(df_pre$t_idx), linetype = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Time",
        y = "Outcome",
        title = "Observed vs Counterfactual (pre & post)"
      ) +
      ggplot2::scale_x_continuous(breaks = cf$t_idx, labels = cf$time)

    return(gg)
  }

  # 以降は従来のタイプ別出力
  td <- tidy_scspill(x)

  if (type == "effect") {
    gg <- ggplot2::ggplot(td$treat, ggplot2::aes(time, mean)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.20) +
      ggplot2::geom_line() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Time (post)",
        y = "Treatment effect",
        title = "Estimated treatment effect (post)"
      )
    return(gg)
  }

  if (type == "spill_top") {
    agg <- stats::aggregate(abs(mean) ~ unit, data = td$spill, FUN = mean)
    top_units <- head(
      agg[order(-agg$`abs(mean)`), "unit"],
      n = min(top_n, nrow(agg))
    )
    df <- td$spill[td$spill$unit %in% top_units, ]
    gg <- ggplot2::ggplot(df, ggplot2::aes(time, mean)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~unit, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Time (post)",
        y = "Spillover effect",
        title = sprintf("Spillover effects (top %d units)", length(top_units))
      )
    return(gg)
  }

  if (type == "weights") {
    ord <- td$weights[order(-abs(td$weights$alpha)), ]
    gg <- ggplot2::ggplot(
      ord,
      ggplot2::aes(x = stats::reorder(unit, abs(alpha)), y = alpha)
    ) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Control unit",
        y = "Synthetic weight (alpha)",
        title = "Estimated synthetic weights"
      )
    return(gg)
  }

  if (type == "rho") {
    df <- data.frame(rho = x$rho_draws)
    gg <- ggplot2::ggplot(df, ggplot2::aes(rho)) +
      ggplot2::geom_histogram(bins = 40) +
      ggplot2::geom_vline(xintercept = mean(df$rho), linetype = 2) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = expression(rho),
        y = "Frequency",
        title = expression(paste("Posterior of ", rho))
      )
    return(gg)
  }

  if (type == "beta") {
    if (is.null(x$sar$beta)) {
      stop("beta draws are not available in `x$sar$beta`.")
    }
    bm <- colMeans(x$sar$beta)
    df <- data.frame(
      name = names(bm) %||% paste0("beta_", seq_along(bm)),
      mean = as.numeric(bm),
      sd = apply(x$sar$beta, 2, sd),
      q025 = apply(x$sar$beta, 2, stats::quantile, 0.025),
      q975 = apply(x$sar$beta, 2, stats::quantile, 0.975)
    )
    gg <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = stats::reorder(name, mean), y = mean)
    ) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = q025, ymax = q975)) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Coefficient",
        y = "Posterior mean (95% CI)",
        title = "Posterior summaries of beta"
      )
    return(gg)
  }

  if (type == "trace") {
    df <- data.frame(iter = seq_along(x$rho_draws), rho = x$rho_draws)
    gg <- ggplot2::ggplot(df, ggplot2::aes(iter, rho)) +
      ggplot2::geom_line() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Iteration",
        y = expression(rho),
        title = expression(paste("Trace of ", rho))
      )
    return(gg)
  }

  stop("Unknown type")
}

# ---------------- S3: autoplot.scspill（ggplot2 流の呼び出し） ----------------

#' @export
#' @method autoplot scspill
autoplot.scspill <- function(object, ...) {
  plot.scspill(object, ...)
}

# ---------------- S3: pp_check.scspill（簡易） ----------------
# ・rho の収束と、post 末期における treated effect の分布近似を提示

#' @export
#' @method pp_check scspill
pp_check.scspill <- function(object, what = c("rho_trace", "treat_last_dist")) {
  what <- match.arg(what)
  if (what == "rho_trace") {
    df <- data.frame(iter = seq_along(object$rho_draws), rho = object$rho_draws)
    gg <- ggplot2::ggplot(df, ggplot2::aes(iter, rho)) +
      ggplot2::geom_line() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Iteration",
        y = expression(rho),
        title = expression(paste("Trace of ", rho))
      )
    return(gg)
  } else {
    mu <- utils::tail(object$effects$treat$mean, 1)
    lo <- utils::tail(object$effects$treat$lo, 1)
    hi <- utils::tail(object$effects$treat$hi, 1)
    sd_ <- (hi - lo) / (2 * 1.96)
    xs <- seq(mu - 4 * sd_, mu + 4 * sd_, length.out = 400)
    df <- data.frame(x = xs, d = stats::dnorm(xs, mean = mu, sd = sd_))
    gg <- ggplot2::ggplot(df, ggplot2::aes(x, d)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Effect (last post period)",
        y = "Density",
        title = "Posterior predictive (normal approx.) vs 0"
      )
    return(gg)
  }
}
