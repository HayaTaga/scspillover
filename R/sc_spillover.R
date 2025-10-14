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