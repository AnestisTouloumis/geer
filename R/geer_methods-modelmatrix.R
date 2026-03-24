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
#' The design matrix is constructed based on \code{terms(object)}, using the data
#' stored in \code{object$data}.
#'
#' For interaction terms, the variable whose levels vary fastest is the first one
#' listed in the formula (not in the term). For example, in
#' \code{~ a + b + b:a}, the interaction term \code{b:a} will vary fastest with
#' respect to \code{a}.
#'
#' By convention, if the response variable also appears on the right-hand side of
#' the formula it is dropped (with a warning). However, interactions involving
#' that term are retained.
#'
#' @return
#' The design matrix for the marginal regression model with the specified formula and
#' data.
#'
#' There is an attribute \code{"assign"}, an integer vector with an entry for
#' each column in the matrix giving the term in the formula which gave rise to
#' the column. Value 0 corresponds to the intercept (if any), and positive
#' values to terms in the order given by the \code{term.labels} attribute of the
#' terms structure corresponding to \code{object}.
#'
#' If there are any factors in terms in the model, there is an attribute
#' \code{"contrasts"}, a named list with an entry for each factor. This
#' specifies the contrasts that would be used in terms in which the factor is
#' coded by contrasts (in some terms dummy coding may be used), either as a
#' character vector naming a function or as a numeric matrix.
#'
#' @examples
#' data("leprosy")
#' fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
#'              family = poisson(link = "log"), id = id, data = leprosy)
#' model.matrix(fit)
#'
#' data("cerebrovascular")
#' fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
#'                     id = id, data = cerebrovascular)
#' model.matrix(fit)
#'
#' @export
model.matrix.geer <-	function(object,...){
  if (!inherits(object, "geer")) {
    stop("'object' must be a 'geer' object", call. = FALSE)
  }
  model.matrix(object = object$terms,
               data = object$data,
               contrasts = object$contrasts)
}
