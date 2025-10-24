#' Horseshoe-Gibbs for synthetic weights alpha
#' @keywords internal
hs_alpha_gibbs <- function(y, X, M = 2000, burn = 1000, verbose = TRUE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  draws <- hs_alpha_gibbs_cpp(y, X, M, burn, verbose)
  colnames(draws) <- paste0("alpha_", seq_len(ncol(X)))
  list(alpha = draws)
}
