.onAttach <- function(lib, pkg) {
  packageStartupMessage("scspill: Bayesian synthetic control with spillovers")
}

#' @useDynLib scspill, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
