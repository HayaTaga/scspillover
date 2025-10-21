#' @keywords internal
scspill_gir_prior_draw <- function(W, N, K, p, priors) {
  scspill_gir_prior_draw_cpp(W, as.integer(N), as.integer(K), as.integer(p),
                             c_beta = priors$c_beta,
                             c_lambda = priors$c_lambda,
                             a0 = priors$a0, b0 = priors$b0)
}

#' @keywords internal
scspill_gir_data_sim <- function(draw, T0, W, w) {
  stopifnot(is.matrix(W), is.numeric(w))
  scspill_gir_data_sim_cpp(draw, as.integer(T0), W, as.numeric(w))
}

#' @keywords internal
scspill_gir_sar_step <- function(theta, sim, W, w, priors) {
  scspill_gir_sar_step_cpp(theta, sim, W, as.numeric(w),
                           c_beta = priors$c_beta, c_lambda = priors$c_lambda,
                           a0 = priors$a0, b0 = priors$b0)
}