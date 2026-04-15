#' @title
#' Construct Design Matrices from a geer Object
#'
#' @aliases model.matrix model.matrix.geer
#' @method model.matrix geer
#'
#' @description
#' Constructs the design matrix for the marginal mean model underlying a fitted
#' \code{geer} object.
#'
#' @inheritParams coef.geer
#'
#' @details
#' The design matrix is reconstructed from \code{object$terms} and
#' \code{object$data}, with factor variables and interactions expanded using
#' the contrasts recorded at fitting time.
#'
#' By convention, if the response variable also appears on the right-hand side
#' of the formula, it is dropped, although interactions involving that term are
#' retained.
#'
#' @return
#' A numeric design matrix for the marginal mean model represented by
#' \code{object}.
#'
#' The returned matrix has an \code{"assign"} attribute: an integer vector with
#' one entry per column giving the index of the formula term that produced that
#' column. A value of \code{0} corresponds to the intercept, if present.
#' Positive values index terms in the order given by
#' \code{attr(terms(object), "term.labels")}.
#'
#' If the model contains factor terms, the returned matrix may also have a
#' \code{"contrasts"} attribute giving the contrasts used for those factors.
#'
#' @seealso \code{\link{geewa}}, \code{\link{geewa_binary}},
#'   \code{\link[stats]{model.matrix}}.
#'
#' @examples
#' data("leprosy", package = "geer")
#' fit <- geewa(
#'   formula = bacilli ~ factor(period) + factor(period):treatment,
#'   family = poisson(link = "log"),
#'   data = leprosy,
#'   id = id
#' )
#' model.matrix(fit)
#'
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id
#' )
#' model.matrix(fit)
#'
#' @export
model.matrix.geer <- function(object, ...) {
  object <- check_geer_object(object)
  model.matrix(
    object = object$terms,
    data = object$data,
    contrasts.arg = object$contrasts
  )
}
