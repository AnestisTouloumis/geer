#' Returns the variance-covariance matrix of the main parameters of a fitted
#' model geer object.
#'
#'
#' Default is to obtain the estimated sandwich (robust) covariance matrix. Setting
#' \code{type = "naive"} obtains the estimated model-based (naive) covariance
#' matrix and \code{type = "bias_corrected"} obtains the estimated bias-corrected
#' covariance matrix.
#'
#' @title Calculate Variance-Covariance Matrix for a Fitted \code{geer} Object.
#'
#' @aliases vcov vcov.geer
#'
#' @method vcov geer
#'
#' @param object a fitted model geer object.
#' @param type character indicating whether the sandwich (robust) covariance
#' matrix (\code{type = "robust"}), the model-based (naive) covariance
#' matrix (\code{type = "naive"}) or the bias-corrected covariance
#' matrix (\code{type = "bias-corrected"}) should be returned.
#' @param ... additional argument(s) for methods.
#'
#' @return A matrix of the estimated covariances between the parameter estimates
#' in the linear predictor of the GEE model. This should have row and column
#' names corresponding to the parameter names given by the coef method.
#'
#'
#' @export

vcov.geer <- function(object, type = "robust", ...) {
  icheck <- pmatch(type,
                   c("robust", "naive", "bias-corrected"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown method for the covariance matrix")
  if (type == "robust") {
    object$robust_covariance
    } else if (type == "bias-corrected") {
      object$bias_corrected_covariance
    } else {
      object$naive_covariance
    }
}
