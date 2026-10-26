#' @title
#' Extract the Dispersion Scale from a geer Object
#'
#' @description
#' Returns the square root of the fitted working dispersion parameter.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param ... not used. Supplying any argument is an error.
#'
#' @details
#' For a Gaussian model, this is the working residual standard deviation.
#' For other families, it is the square root of the dispersion multiplying
#' the variance function. For \code{geewa_binary()} fits, it is 1.
#' The method also supplies the Gaussian scale used by \pkg{emmeans}.
#'
#' @return
#' A positive numeric scalar.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{emmeans-support}}.
#'
#' @examples
#' data("rinse", package = "geer")
#' fit <- geewa(
#'   formula = score ~ treatment + baseline + time,
#'   data = rinse,
#'   id = id
#' )
#' sigma(fit)
#'
#' @export
sigma.geer <- function(object, ...) {
  check_unused_dots(list(...), "sigma.geer")
  object <- check_geer_object(object)
  phi <- check_single_numeric(object$phi, "object$phi")
  if (phi <= 0) {
    stop("'object$phi' must be positive", call. = FALSE)
  }
  sqrt(phi)
}
