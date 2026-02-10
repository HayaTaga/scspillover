.scspill_required_inputs <- c(
  "Y0_pre",
  "Y0_post",
  "Yc_pre",
  "Yc_post",
  "times_pre",
  "times_post",
  "units"
)

#' @keywords internal
validate_scspill <- function(x) {
  stopifnot(is.list(x), inherits(x, "scspill"))

  miss_top <- setdiff(.scspill_required_inputs, names(x$inputs))
  if (length(miss_top)) {
    stop(
      "scspill object is missing top-level inputs: ",
      paste(miss_top, collapse = ", ")
    )
  }
  if (is.null(x$inputs$units$control)) {
    stop("scspill object is missing inputs$units$control")
  }

  T0 <- nrow(as.matrix(x$inputs$Yc_pre))
  T1 <- nrow(as.matrix(x$inputs$Yc_post))
  Np <- ncol(as.matrix(x$inputs$Yc_pre))
  N1 <- ncol(as.matrix(x$inputs$Yc_post))

  if (length(x$inputs$Y0_pre) != T0) {
    stop("Y0_pre length mismatch.")
  }
  if (length(x$inputs$Y0_post) != T1) {
    stop("Y0_post length mismatch.")
  }
  if (length(as.vector(x$inputs$times_pre)) != T0) {
    stop("times_pre length mismatch.")
  }
  if (length(as.vector(x$inputs$times_post)) != T1) {
    stop("times_post length mismatch.")
  }
  if (length(as.character(x$inputs$units$control)) != Np || Np != N1) {
    stop("units$control length / Yc_pre/post ncol mismatch.")
  }

  invisible(TRUE)
}

#' @keywords internal
new_scspill <- function(
  alpha_draws,
  rho_draws,
  alpha_hat,
  rho_hat,
  effects,
  inputs,
  sar = NULL
) {
  inputs$Yc_pre <- as.matrix(inputs$Yc_pre)
  storage.mode(inputs$Yc_pre) <- "double"
  inputs$Yc_post <- as.matrix(inputs$Yc_post)
  storage.mode(inputs$Yc_post) <- "double"
  inputs$Y0_pre <- as.numeric(inputs$Y0_pre)
  inputs$Y0_post <- as.numeric(inputs$Y0_post)

  if (is.null(inputs$times_pre)) {
    T0 <- nrow(inputs$Yc_pre)
    inputs$times_pre <- if (!is.null(rownames(inputs$Yc_pre))) {
      rownames(inputs$Yc_pre)
    } else {
      seq_len(T0)
    }
  }
  if (is.null(inputs$times_post)) {
    T1 <- nrow(inputs$Yc_post)
    inputs$times_post <- if (!is.null(rownames(inputs$Yc_post))) {
      rownames(inputs$Yc_post)
    } else {
      seq_len(T1)
    }
  }
  inputs$times_pre <- as.vector(inputs$times_pre)
  inputs$times_post <- as.vector(inputs$times_post)

  if (is.null(inputs$units) || is.null(inputs$units$control)) {
    u <- colnames(inputs$Yc_post)
    if (is.null(u)) {
      u <- colnames(inputs$Yc_pre)
    }
    if (is.null(inputs$units)) {
      inputs$units <- list()
    }
    inputs$units$control <- if (!is.null(u)) {
      as.character(u)
    } else {
      paste0("unit_", seq_len(ncol(inputs$Yc_pre)))
    }
  }

  if (!is.null(inputs$w)) {
    inputs$w <- as.matrix(inputs$w)
  }
  if (!is.null(inputs$W)) {
    inputs$W <- as.matrix(inputs$W)
  }

  obj <- structure(
    list(
      alpha_draws = as.matrix(alpha_draws),
      rho_draws = as.numeric(rho_draws),
      alpha_hat = as.numeric(alpha_hat),
      rho_hat = as.numeric(rho_hat),
      effects = effects,
      sar = sar,
      inputs = inputs
    ),
    class = "scspill"
  )

  validate_scspill(obj)
  obj
}
