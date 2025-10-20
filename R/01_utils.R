#' Prepare panel & covariates into matrices (controls only)
#' @keywords internal
scspill_prep_X <- function(
  data,
  treated_unit,
  T0,
  X = NULL,
  y_col = "y",
  unit_col = "unit",
  time_col = "time"
) {
  stopifnot(all(c(unit_col, time_col, y_col) %in% names(data)))

  d <- data
  # 内部表記統一
  names(d)[names(d) == unit_col] <- "unit"
  names(d)[names(d) == time_col] <- "time"
  if (y_col != "y") {
    d$y <- d[[y_col]]
  }

  d <- d[order(d$time, d$unit), ]
  times <- sort(unique(d$time))
  T <- length(times)
  units <- sort(unique(d$unit))
  idx0 <- which(units == treated_unit)
  controls <- units[-idx0]
  N <- length(controls)

  Ywide <- d %>%
    pivot_wider(
      id_cols = time,
      names_from = unit,
      values_from = y
    ) %>%
    arrange(time)
  Ywide <- Ywide[order(Ywide$time), ]
  Y <- as.matrix(Ywide[, -1])
  colnames(Y) <- sub("^y\\.", "", colnames(Y))

  Y0 <- Y[, treated_unit, drop = TRUE]
  Yc <- Y[, controls, drop = FALSE]

  X_3d <- NULL
  if (!is.null(X)) {
    if (is.character(X)) {
      # y_col が含まれていれば自動除外
      X <- setdiff(X, y_col)
      stopifnot(all(X %in% names(d)))
      Xmat <- as.matrix(d[, X, drop = FALSE])
    } else if (is.matrix(X) || is.data.frame(X)) {
      Xmat <- as.matrix(X)
      stopifnot(nrow(Xmat) == nrow(d))
    } else {
      stop(
        "X must be NULL, character colnames in data, or a matrix/data.frame."
      )
    }
    K <- ncol(Xmat)
    X_3d <- array(NA_real_, c(T, N, K))
    for (ti in seq_len(T)) {
      rows_t <- which(d$time == times[ti])
      df_t <- d[rows_t, , drop = FALSE]
      Xm_t <- Xmat[rows_t, , drop = FALSE]
      ord <- match(controls, df_t$unit)
      X_3d[ti, , ] <- Xm_t[ord, , drop = FALSE]
    }
  }

  list(
    Y0_pre = Y0[1:T0],
    Yc_pre = Yc[1:T0, , drop = FALSE],
    Y0_post = Y0[(T0 + 1):T],
    Yc_post = Yc[(T0 + 1):T, , drop = FALSE],
    Xc_pre = if (is.null(X_3d)) NULL else X_3d[1:T0, , , drop = FALSE],
    Xc_post = if (is.null(X_3d)) NULL else X_3d[(T0 + 1):T, , , drop = FALSE],
    times = times,
    units = list(treated = treated_unit, controls = controls)
  )
}

#' Posterior effects via identification formulas (5)(6)
#' @keywords internal
posterior_effects <- function(
  Y0_post,
  Yc_post,
  alpha_draws,
  rho_draws,
  w,
  W,
  cred = 0.95
) {
  T1 <- length(Y0_post)
  N <- ncol(Yc_post)
  M <- nrow(alpha_draws)
  Aeff <- array(NA_real_, dim = c(T1, M))
  Spill <- array(NA_real_, dim = c(T1, N, M))
  IN <- diag(N)

  for (m in 1:M) {
    a <- as.matrix(alpha_draws[m, ])
    r <- rho_draws[m]
    # Ainv <- solve(IN - r * (w %*% t(a) + W)) # (IN - r w a' - r W)^{-1}
    A <- inverse_check(IN, r, w, a, W)      # 固有値半径をもとに r を安全化し、A を組み立て
    Ainv <- robust_solve(A) 
    B <- (IN - r * W)
    for (t in 1:T1) {
      yc <- Yc_post[t, ]
      y0 <- Y0_post[t]
      tmp <- Ainv %*% (B %*% yc - r * w * y0)
      Aeff[t, m] <- y0 - as.numeric(crossprod(a, tmp)) # (5)
      Spill[t, , m] <- yc - as.vector(tmp) # (6)
    }
  }
  treat_mean <- rowMeans(Aeff)
  treat_q <- apply(Aeff, 1, quantile, probs = c(0.025, 0.975))
  spill_mean <- apply(Spill, c(1, 2), mean)
  spill_lo <- apply(Spill, c(1, 2), quantile, probs = 0.025)
  spill_hi <- apply(Spill, c(1, 2), quantile, probs = 0.975)

  list(
    treat = list(mean = treat_mean, lo = treat_q[1, ], hi = treat_q[2, ]),
    spill = list(mean = spill_mean, lo = spill_lo, hi = spill_hi)
  )
}

inverse_check <- function(IN, r, w, a, W, eps_spec = 1e-3) {
  N <- nrow(IN)
  stopifnot(ncol(IN) == N, length(w) == N, length(a) == ncol(W), nrow(W) == N)

  B <- w %*% t(a) + W

  ev  <- eigen(B, only.values = TRUE)$values
  rho <- max(Mod(ev))
  if (is.finite(rho) && rho > 0) {
    r_max <- (1 - eps_spec) / rho
    if (abs(r) >= r_max) r <- sign(r) * r_max
  }

  A <- IN - r * B
  A
}

robust_solve <- function(A, b = NULL, ridge0 = 1e-12, max_tries = 6, qr_tol = 1e-12) {
  N <- nrow(A); I <- diag(N)
  lam <- 0
  for (k in 0:max_tries) {
    Areg <- if (lam == 0) A else A + lam * I
    # rcond/kappa の評価（失敗時は次へ）
    ok <- tryCatch({
      rc <- 1 / kappa(Areg, exact = FALSE)
      is.finite(rc) && rc > 1e-12
    }, error = function(e) FALSE)

    if (ok) {
      # 右辺があれば Ax=b を解き、無ければ擬似逆行列的に I を右辺に
      if (is.null(b)) {
        # 逆行列が本当に必要なら列ごとに解く方が安定
        return(qr.solve(Areg, diag(N), tol = qr_tol))
      } else {
        return(qr.solve(Areg, b, tol = qr_tol))
      }
    }
    lam <- if (lam == 0) ridge0 else lam * 10
  }
  # 最後の手段（MASS::ginv）。理論的に特異で許容できる場合にのみ。
  if (!requireNamespace("MASS", quietly = TRUE)) stop("MASS not installed for ginv fallback.")
  if (is.null(b)) return(MASS::ginv(A))
  MASS::ginv(A) %*% b
}