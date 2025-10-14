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
              x$rho_mean, x$rho_ci[1], x$rho_ci[2]))
  cat("non-negligible alpha indices:",
      if (length(x$alpha_nonzero)) paste(x$alpha_nonzero, collapse = ", ") else "(none)", "\n")
  cat(sprintf("average post-treatment effect (treated unit) = %.3f\n", x$treat_avg))
}

#' @export
plot.scspill <- function(x, type = c("treated", "spillover"), unit_index = NULL) {
  type <- match.arg(type)
  T1 <- length(x$effects$treat$mean)
  tt <- seq_len(T1)
  if (type == "treated") {
    df <- data.frame(t = tt, m = x$effects$treat$mean,
                     lo = x$effects$treat$lo, hi = x$effects$treat$hi)
    ggplot2::ggplot(df, ggplot2::aes(t, m)) + ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Post-treatment time", y = "Treatment effect")
  } else {
    if (is.null(unit_index)) stop("unit_index is required for type='spillover'")
    df <- data.frame(t = tt,
                     m = x$effects$spill$mean[, unit_index],
                     lo = x$effects$spill$lo[, unit_index],
                     hi = x$effects$spill$hi[, unit_index])
    ggplot2::ggplot(df, ggplot2::aes(t, m)) + ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
      ggplot2::theme_minimal() + ggplot2::labs(x = "Post-treatment time", y = "Spillover effect")
  }
}