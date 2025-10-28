# ---------------- 共通ユーティリティ ----------------

`%||%` <- function(x, y) if (!is.null(x)) x else y

.as_num <- function(x) {
  if (is.null(x)) {
    return(numeric(0))
  }
  as.numeric(x)
}
.as_mat <- function(x) {
  if (is.null(x)) {
    return(matrix(numeric(0), 0, 0))
  }
  as.matrix(x)
}

.safe_seq <- function(n, start = 1L) {
  n <- as.integer(n)
  start <- as.integer(start)
  if (is.na(n) || n <= 0L) integer(0) else seq.int(from = start, length.out = n)
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

# ----------------------- times 取得（time_col 対応） -----------------------
.get_times_pre <- function(fit, time_col = NULL) {
  # data_pre の列を最優先
  if (
    !is.null(time_col) &&
      !is.null(fit$inputs$data_pre) &&
      is.data.frame(fit$inputs$data_pre) &&
      time_col %in% names(fit$inputs$data_pre)
  ) {
    return(as.vector(fit$inputs$data_pre[[time_col]]))
  }
  # 既存フィールド
  if (!is.null(fit$inputs$times_pre)) {
    return(as.vector(fit$inputs$times_pre))
  }

  rn <- rownames(fit$inputs$Yc_pre)
  if (!is.null(rn)) {
    if (suppressWarnings(all(!is.na(as.numeric(rn))))) {
      return(as.numeric(rn))
    }
    return(rn)
  }
  seq_len(nrow(fit$inputs$Yc_pre))
}

.get_times_post <- function(fit, time_col = NULL) {
  if (
    !is.null(time_col) &&
      !is.null(fit$inputs$data_post) &&
      is.data.frame(fit$inputs$data_post) &&
      time_col %in% names(fit$inputs$data_post)
  ) {
    return(as.vector(fit$inputs$data_post[[time_col]]))
  }
  if (!is.null(fit$inputs$times_post)) {
    return(as.vector(fit$inputs$times_post))
  }
  rn <- rownames(fit$inputs$Yc_post)
  if (!is.null(rn)) {
    if (suppressWarnings(all(!is.na(as.numeric(rn))))) {
      return(as.numeric(rn))
    }
    return(rn)
  }
  seq_len(nrow(fit$inputs$Yc_post))
}

.standardize_spill <- function(spill_raw, times_post, unit_names = NULL) {
  if (is.null(spill_raw)) {
    return(data.frame())
  } # なし

  # list(mean=matrix, ...) も扱う
  if (is.list(spill_raw) && !"data.frame" %in% class(spill_raw)) {
    if (!is.null(spill_raw$mean)) spill_raw <- spill_raw$mean
  }

  # case 1: すでにロング（unit/time/mean を持つ）
  if (
    is.data.frame(spill_raw) &&
      all(c("unit", "time", "mean") %in% names(spill_raw))
  ) {
    df <- spill_raw[, c("unit", "time", "mean")]
    df$unit <- as.character(df$unit)
    df$time <- suppressWarnings(as.numeric(df$time))
    df$mean <- suppressWarnings(as.numeric(df$mean))
    df <- df[
      !is.na(df$unit) & !is.na(df$time) & !is.na(df$mean),
      ,
      drop = FALSE
    ]
    return(df)
  }

  # case 2: 行列・ワイド（列＝ユニット）
  if (
    is.matrix(spill_raw) ||
      (is.data.frame(spill_raw) && !("unit" %in% names(spill_raw)))
  ) {
    M <- as.matrix(spill_raw)
    T1 <- nrow(M)
    U <- ncol(M)
    if (U <= 1) {
      return(data.frame())
    } # ユニット別がない場合は空で返す
    # time
    if (length(times_post) != T1) {
      times_post <- .safe_seq(T1, 1L)
    }
    # unit 名
    coln <- colnames(M)
    if (is.null(coln) || any(coln == "")) {
      if (!is.null(unit_names)) {
        coln <- unit_names[seq_len(U)]
      } else {
        coln <- paste0("unit_", seq_len(U))
      }
    }
    # ロングへ（列優先ベクトル化に合わせて rep の向きに注意）
    df <- data.frame(
      time = rep(times_post, times = U),
      unit = rep(coln, each = T1),
      mean = as.numeric(M)
    )
    df$unit <- as.character(df$unit)
    df$time <- suppressWarnings(as.numeric(df$time))
    df$mean <- suppressWarnings(as.numeric(df$mean))
    df <- df[!is.na(df$mean), , drop = FALSE]
    return(df)
  }

  # それ以外は扱えない
  data.frame()
}

# -------------- counterfactual（事後平均/区間を計算） ----------------
#   Ainv = (I_N - r (w a' + W))^{-1}
#   tmp  = Ainv { (I_N - r W) y_c - r w y0 }
#   effect = y0 - a' tmp  =>  counterfactual = a' tmp
#
# 返り値: data.frame(time, t_idx, period, y_obs, y_cf_mean, y_cf_lo, y_cf_hi)
scspill_counterfactual <- function(fit, cred = 0.95, time_col = NULL) {
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

  # 時間軸（指定があれば data_* の列を使う）
  times_pre <- .get_times_pre(fit, time_col)
  if (length(times_pre) != T0) {
    times_pre <- .safe_seq(T0, 1L)
  }
  times_post <- .get_times_post(fit, time_col)
  if (length(times_post) != T1) {
    times_post <- .safe_seq(T1, if (T0 > 0) times_pre[T0] + 1L else 1L)
  }

  # pre
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

  # post
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

  df_pre <- data.frame(
    time = times_pre,
    t_idx = seq_len(T0),
    period = "pre",
    y_obs = Y0_pre,
    y_cf_mean = rowMeans(ycf_pre_draws),
    y_cf_lo = apply(ycf_pre_draws, 1, stats::quantile, probs = lo_q),
    y_cf_hi = apply(ycf_pre_draws, 1, stats::quantile, probs = hi_q)
  )
  df_post <- data.frame(
    time = times_post,
    t_idx = T0 + seq_len(T1),
    period = "post",
    y_obs = Y0_post,
    y_cf_mean = rowMeans(ycf_post_draws),
    y_cf_lo = apply(ycf_post_draws, 1, stats::quantile, probs = lo_q),
    y_cf_hi = apply(ycf_post_draws, 1, stats::quantile, probs = hi_q)
  )
  rbind(df_pre, df_post)
}

# ---------------- tidy 化（効果・スピル・重み等を取り出し） ----------------
# 既存の fit$effects / fit$weights があればそれを尊重。無ければ最小限を再構成。
tidy_scspill <- function(fit, time_col = NULL) {
  y0_pre <- as.numeric(fit$inputs$Y0_pre)
  y0_post <- as.numeric(fit$inputs$Y0_post)
  Yc_pre <- as.matrix(fit$inputs$Yc_pre)
  Yc_post <- as.matrix(fit$inputs$Yc_post)

  T0 <- length(y0_pre)
  T1 <- length(y0_post)
  T <- T0 + T1

  # post の time（spill でも使う）
  times_post <- .get_times_post(fit, time_col)
  if (length(times_post) != T1) {
    times_post <- .safe_seq(T1, if (T0 > 0) T0 + 1L else 1L)
  }

  # コントロール名（spill の列名補完に使う）
  unit_names <- NULL
  if (!is.null(fit$inputs$units) && is.list(fit$inputs$units)) {
    unit_names <- as.character(fit$inputs$units$control)
  } else {
    unit_names <- colnames(Yc_post) %||% colnames(Yc_pre)
  }

  weights_df <- data.frame()
  if (!is.null(fit$alpha_draws)) {
    A <- as.matrix(fit$alpha_draws) # M x N
    if (ncol(A) > 0) {
      alpha_hat <- colMeans(A, na.rm = TRUE)
      units <- .get_units_control(fit, length(alpha_hat))
      # ここで必ず numeric に落とす
      weights_df <- data.frame(
        unit = as.character(units),
        alpha = as.numeric(alpha_hat),
        stringsAsFactors = FALSE
      )
    }
  }

  # spill をロング化
  spill_df <- data.frame()
  if (!is.null(fit$effects) && !is.null(fit$effects$spill)) {
    spill_df <- .standardize_spill(fit$effects$spill, times_post, unit_names)
  }

  list(
    T0 = T0,
    T1 = T1,
    T = T,
    times_post = times_post,
    y0_all = c(y0_pre, y0_post),
    Yc_all = rbind(Yc_pre, Yc_post),
    spill = spill_df,
    weights = weights_df
  )
}

# ---------------- S3: plot.scspill ----------------

#' scspill の可視化
#'
#' @param x scspill オブジェクト
#' @param type "full"（既定: 実測vs反事実, 前後一括）, "effect"（postの効果推移）,
#'             "spill_top", "weights", "rho", "beta", "trace"
#' @param top_n spill_top の表示ユニット数
#' @param cred  信用水準（反事実リボン）
#' @param time_col data_pre / data_post にある横軸列名（例 "year"）
#' @param ...   予備
#' @export
plot.scspill <- function(
  x,
  type = c("full", "effect", "spill_top", "weights", "rho", "beta", "trace"),
  top_n = 8,
  cred = 0.95,
  time_col = NULL,
  ...
) {
  type <- match.arg(type)

  # -------- full: Observed vs Counterfactual --------
  if (type == "full") {
    cf <- scspill_counterfactual(x, cred = cred, time_col = time_col)

    # 連番 index（pre→post で単調増加。線は idx で結び、ラベルは time）
    cf$.idx <- seq_len(nrow(cf))
    T0 <- sum(cf$period == "pre")

    df_pre <- subset(cf, period == "pre")
    df_post <- subset(cf, period == "post")

    # CF 線を介入年（pre 最終点）から連続させる
    if (nrow(df_pre) > 0 && nrow(df_post) > 0) {
      pre_last <- df_pre[nrow(df_pre), ]
      pre_last$period <- "post"
      df_post_line <- rbind(pre_last, df_post)
    } else {
      df_post_line <- df_post
    }

    # 線データ（描画は .idx、ラベルは time）
    obs_line <- data.frame(idx = cf$.idx, y = cf$y_obs, series = "Observed")
    cf_line_pre <- data.frame(
      idx = df_pre$.idx,
      y = df_pre$y_cf_mean,
      series = "Counterfactual"
    )
    cf_line_post <- data.frame(
      idx = df_post_line$.idx,
      y = df_post_line$y_cf_mean,
      series = "Counterfactual"
    )

    # リボン用
    rib_pre <- transform(df_pre, idx = .idx)
    rib_post <- transform(df_post, idx = .idx)

    # x 軸（ラベルは time）
    x_breaks <- cf$.idx
    x_labels <- cf$time

    gg <- ggplot2::ggplot() +
      # 観測
      ggplot2::geom_line(
        data = obs_line,
        ggplot2::aes(idx, y, linetype = series, color = series)
      ) +
      # CF リボン（pre / post）
      ggplot2::geom_ribbon(
        data = rib_pre,
        ggplot2::aes(
          x = idx,
          ymin = y_cf_lo,
          ymax = y_cf_hi,
          fill = "Counterfactual 95% CI"
        ),
        alpha = 0.12,
        show.legend = TRUE
      ) +
      ggplot2::geom_ribbon(
        data = rib_post,
        ggplot2::aes(
          x = idx,
          ymin = y_cf_lo,
          ymax = y_cf_hi,
          fill = "Counterfactual 95% CI"
        ),
        alpha = 0.18,
        show.legend = TRUE
      ) +
      # CF 線（介入年から連続）
      ggplot2::geom_line(
        data = cf_line_pre,
        ggplot2::aes(idx, y, linetype = series, color = series)
      ) +
      ggplot2::geom_line(
        data = cf_line_post,
        ggplot2::aes(idx, y, linetype = series, color = series)
      ) +
      # 介入境界（pre の最後）
      ggplot2::geom_vline(xintercept = T0, linetype = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Time",
        y = "Outcome",
        title = "Observed vs Counterfactual (pre & post)"
      ) +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
      # 凡例：上・横並び（線種で区別、色は見た目だけ・凡例非表示）
      ggplot2::scale_linetype_manual(
        values = c("Observed" = "solid", "Counterfactual" = "22"),
        name = NULL
      ) +
      ggplot2::scale_color_manual(
        values = c("Observed" = "black", "Counterfactual" = "black"),
        guide = "none",
        name = NULL
      ) +
      ggplot2::scale_fill_manual(
        values = c("Counterfactual 95% CI" = "grey70"),
        name = NULL
      ) +
      ggplot2::theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::guides(
        linetype = ggplot2::guide_legend(order = 1, nrow = 1),
        fill = ggplot2::guide_legend(order = 2, nrow = 1)
      )

    return(gg)
  }

  # 以降は tidy 済みオブジェクトを利用
  td <- tidy_scspill(x, time_col = time_col)

  # -------- effect: post の処置効果 --------
  if (type == "effect") {
    
    # CFを計算（"full" と同じロジック）
    cf <- scspill_counterfactual(x, cred = cred, time_col = time_col)
    
    # 介入後 (post) のみ抽出
    df_post <- subset(cf, period == "post")
    
    # 処置効果 (Effect = Observed - Counterfactual) を計算
    # 信用区間は (Obs - CF_hi, Obs - CF_lo) となる
    df <- data.frame(
        time = df_post$time,
        mean = df_post$y_obs - df_post$y_cf_mean,
        lo = df_post$y_obs - df_post$y_cf_hi, # 観測値 - CFの上限 = 効果の下限
        hi = df_post$y_obs - df_post$y_cf_lo  # 観測値 - CFの下限 = 効果の上限
    )
    
    # エラー回避: 介入後データがない場合は空のプロットを返す
    if (nrow(df) == 0) {
        warning("`type = 'effect'` が呼び出されましたが、介入後のデータ (post) が見つかりません。")
        return(ggplot2::ggplot() + ggplot2::theme_void() + 
               ggplot2::ggtitle("No post-treatment data found."))
    }

    # 連続 index（ラベルは time）
    df$.idx <- seq_len(nrow(df))
    
    gg <- ggplot2::ggplot(df, ggplot2::aes(.idx, mean)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.20) +
      # ゼロラインを追加
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_line() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Time (post)",
        y = "Treatment effect (Observed - Counterfactual)", # Y軸ラベルを明確化
        title = "Estimated treatment effect (post)"
      ) +
      ggplot2::scale_x_continuous(breaks = df$.idx, labels = df$time)
    return(gg)
  }
  # -------- spill_top: 上位ユニットのスピル --------
  if (type == "spill_top") {
    df <- td$spill 

    if (!nrow(df)) {
      stop(
        "`effects$spill` にユニット別のスピル系列が見つかりません（単一系列か未保存の可能性）。\n",
        "各ユニット列を持つ行列/データフレーム、または (unit,time,mean) のロング形式を入れてください。"
      )
    }

    df$unit <- as.character(df$unit)
    df$time <- suppressWarnings(as.numeric(df$time))
    df$mean <- suppressWarnings(as.numeric(df$mean))
    df <- df[
      !is.na(df$unit) & !is.na(df$time) & !is.na(df$mean),
      ,
      drop = FALSE
    ]

    agg <- stats::aggregate(
      x = list(abs_mean = abs(df$mean)),
      by = list(unit = df$unit),
      FUN = mean,
      na.rm = TRUE
    )
    agg <- agg[order(-agg$abs_mean), , drop = FALSE]
    top_units <- head(agg$unit, n = min(top_n, nrow(agg)))

    df <- df[df$unit %in% top_units, , drop = FALSE]
    df <- df[order(df$unit, df$time), , drop = FALSE]

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

  # -------- weights --------
  if (type == "weights") {
    td <- tidy_scspill(x, time_col = time_col)
    df <- td$weights

    if (!nrow(df)) {
      stop(
        "weights（alpha の事後要約）が見つかりません。`fit$alpha_draws` が空でないか確認してください。"
      )
    }

    # 念のため型を保証（ここで numeric に）
    df$unit <- as.character(df$unit)
    df$alpha <- suppressWarnings(as.numeric(df$alpha))
    df <- df[!is.na(df$alpha), , drop = FALSE]

    ord <- df[order(-abs(df$alpha)), , drop = FALSE]

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
        y = "Synthetic weight (alpha, posterior mean)",
        title = "Estimated synthetic weights"
      )
    return(gg)
  }

  # -------- rho: 事後分布のヒスト --------
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

  # -------- beta: 係数の要約 --------
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

  # -------- trace: rho のトレース --------
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
    # treat の最後の点を正規近似
    td <- tidy_scspill(object)
    if (nrow(td$treat) == 0) {
      stop("No treatment-effect summary found.")
    }
    mu <- utils::tail(td$treat$mean, 1)
    lo <- utils::tail(td$treat$lo, 1)
    hi <- utils::tail(td$treat$hi, 1)
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
