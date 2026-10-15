#' @title
#' Construct Design Matrices from a geer Object
#'
#' @method model.matrix geer
#'
#' @description
#' Extracts the model matrix of the marginal mean model of a fitted \code{geer}
#' object.
#'
#' @inheritParams coef.geer
#'
#' @details
#' The matrix is the one stored when the model was fitted, so it agrees with
#' the coefficients, contrasts, and observations used in the fit, and its rows
#' are in the internal order of the fitted model.
#'
#' @return
#' The model matrix, with attributes \code{"assign"}, which maps the columns
#' to the terms (0 for the intercept), and, for models with factors,
#' \code{"contrasts"}.
#'
#' @seealso \code{\link{geewa}}, \code{\link{geewa_binary}},
#'   \code{\link[stats]{model.matrix}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' head(model.matrix(fit))
#'
#' @export
model.matrix.geer <- function(object, ...) {
  object <- check_geer_object(object)
  object$x
}
