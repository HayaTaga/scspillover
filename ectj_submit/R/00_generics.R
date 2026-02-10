#' @importFrom ggplot2 autoplot
NULL

#' Posterior predictive check generic
#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}

#' Diagnostics generic
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}
