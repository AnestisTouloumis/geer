#' @title
#' Construct Design Matrices from a \code{geer} Object
#'
#' @aliases model.matrix
#' @method model.matrix geer
#'
#' @description
#' Construct the design (model) matrix for the marginal mean model underlying a
#' fitted \code{geer} object. Factor variables are expanded according to the
#' contrasts used when fitting the model, and interaction terms are expanded
#' accordingly.
#'
#' @inheritParams coef.geer
#'
#' @details
#' The design matrix is reconstructed from the terms object stored in
#' \code{object$terms} and the data stored in \code{object$data}. Factor
#' variables are expanded according to the contrasts used when fitting the
#' model, and interaction terms are expanded accordingly.
#'
#' By convention, if the response variable also appears on the right-hand side
#' of the formula it is dropped, although interactions involving that term are
#' retained.
#'
#' @return
#' A numeric design matrix for the marginal mean model represented by
#' \code{object}.
#'
#' There is an attribute \code{"assign"}, an integer vector with an entry for
#' each column in the matrix giving the term in the formula which gave rise to
#' the column. Value 0 corresponds to the intercept, if present, and positive
#' values correspond to terms in the order given by the \code{term.labels}
#' attribute of \code{terms(object)}.
#'
#' If the model contains factor terms, the returned matrix may also have a
#' \code{"contrasts"} attribute giving the contrasts used for those factors.
#'
#' @examples
#' data("leprosy", package = "geer")
#' fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
#'              family = poisson(link = "log"), id = id, data = leprosy)
#' model.matrix(fit)
#'
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
#'                     id = id, data = cerebrovascular)
#' model.matrix(fit)
#'
#' @export
model.matrix.geer <-	function(object,...){
  object <- check_geer_object(object)
  ans <- model.matrix(
    object = object$terms,
    data = object$data,
    contrasts.arg = object$contrasts
    )
  ans
}
