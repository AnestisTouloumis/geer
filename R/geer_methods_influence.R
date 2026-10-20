#' @title
#' Leverage and Deletion Diagnostics for a geer Object
#'
#' @rdname influence.geer
#' @method hatvalues geer
#'
#' @description
#' Computes leverages, the change in the estimated regression coefficients
#' when a cluster or an observation is deleted, and Cook's distances for a
#' fitted \code{geer} object.
#'
#' @param model a fitted model object of class \code{"geer"}.
#' @param level character string specifying whether the diagnostics refer to
#'   clusters (\code{"clusters"}) or to single observations
#'   (\code{"observations"}). Defaults to \code{"clusters"}.
#' @param method character string specifying how the deletion estimates are
#'   obtained: the one-step approximation of Preisser and Qaqish (1996)
#'   (\code{"one-step"}) or a refit of the model without each cluster
#'   (\code{"exact"}). See Details. Defaults to \code{"one-step"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used to standardize the deletion estimates: one of the estimators
#'   described in \code{\link{vcov.geer}}. Defaults to \code{"naive"}.
#' @param parm a specification of which coefficients enter Cook's distance,
#'   either a vector of numbers or a vector of names. If missing, all
#'   coefficients are used.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @details
#' Let \eqn{D_i} be the derivative of the mean vector of cluster \eqn{i} with
#' respect to \eqn{\beta}, \eqn{V_i} its working covariance matrix, including
#' the dispersion, and
#' \eqn{\Sigma = (\sum_{i} D_i^\top V_i^{-1} D_i)^{-1}}{Sigma =
#' (sum_i D_i' V_i^(-1) D_i)^(-1)} the naive covariance matrix. The leverage
#' matrix of cluster \eqn{i} is
#' \deqn{H_i = D_i \Sigma D_i^\top V_i^{-1}.}{H_i = D_i Sigma D_i' V_i^(-1).}
#' \code{hatvalues()} returns the diagonal elements of \eqn{H_i} for
#' \code{level = "observations"} and their sum, the trace of \eqn{H_i}, for
#' \code{level = "clusters"}. The leverages sum to the number of coefficients.
#'
#' \code{dfbeta()} returns \eqn{\hat\beta - \hat\beta_{(i)}}{beta - beta_(i)},
#' where \eqn{\hat\beta_{(i)}}{beta_(i)} is the estimate without cluster or
#' observation \eqn{i}. With \code{method = "one-step"}, it is a single
#' Fisher scoring step from \eqn{\hat\beta}{beta} with the association and
#' dispersion parameters held at their full-data values (Preisser and Qaqish,
#' 1996). For a cluster this is
#' \deqn{\Sigma D_i^\top V_i^{-1} (I - H_i)^{-1} (y_i - \hat\mu_i),}{Sigma
#' D_i' V_i^(-1) (I - H_i)^(-1) (y_i - mu_i),}
#' as in equation (4.66) of Hardin and Hilbe (2013). For observation \eqn{j}
#' of cluster \eqn{i}, the step uses the residual, the derivative and the
#' variance of \eqn{y_{ij}}{y_ij} conditional on the remaining observations of
#' the cluster under the working covariance matrix. This differs from equation
#' (4.67) of Hardin and Hilbe (2013), which ignores the correlation with the
#' rest of the cluster and agrees with it only under working independence.
#' Both one-step approximations are exact for a Gaussian response with the
#' identity link.
#'
#' With \code{method = "exact"}, every cluster is deleted in turn and the model
#' refitted with the association parameters held at their full-data values,
#' as for \code{cov_type = "jackknife"} in \code{\link{vcov.geer}}. The
#' refits are computed once and cached on the fit. It is available only for
#' \code{level = "clusters"}. The one-step approximation is derived from the
#' ordinary generalized estimating equations, so it is available only for
#' fits with \code{method = "gee"}; for the other estimation methods, use
#' \code{method = "exact"}.
#'
#' When deleting a cluster or an observation leaves coefficients that cannot be
#' estimated, for example the only cluster at a level of a factor, the
#' leverage matrix \eqn{H_i} has an eigenvalue of 1 and the deletion estimate
#' is undefined. Its values are then \code{NA}, with a warning naming the
#' clusters or observations, and no refit is attempted for them with
#' \code{method = "exact"}.
#'
#' \code{dfbetas()} divides each column of \code{dfbeta()} by the standard
#' error of the coefficient from \code{vcov(model, cov_type = cov_type)}.
#' \code{cooks.distance()} returns
#' \deqn{\frac{1}{q} (\hat\beta - \hat\beta_{(i)})^\top \hat V^{-1}
#' (\hat\beta - \hat\beta_{(i)}),}{(1/q) (beta - beta_(i))' V^(-1) (beta -
#' beta_(i)),}
#' where the coefficients are restricted to those in \code{parm}, \eqn{q} is
#' their number and \eqn{\hat V}{V} their covariance matrix from
#' \code{vcov(model, cov_type = cov_type)}. With the default
#' \code{cov_type = "naive"} and all coefficients, this is equation (4.68) of
#' Hardin and Hilbe (2013) for clusters, and the corresponding measure of
#' Preisser and Qaqish (1996) for observations. \pkg{glmtoolbox} uses the
#' robust covariance matrix by default.
#'
#' @section Observation order:
#' For \code{level = "observations"}, the values follow the internal order of
#' the fitted model: the complete cases, sorted by cluster and, within
#' cluster, by \code{repeated}. They are named \code{"id(repeated)"} by the
#' cluster identifier and the value of \code{repeated}, or the position within
#' the cluster if \code{repeated} was not supplied.
#'
#' @return
#' For \code{hatvalues()} and \code{cooks.distance()}, a numeric vector with
#' one value per cluster, named by the cluster identifier, or one value per
#' observation, in the order and with the names given under Observation
#' order. For
#' \code{dfbeta()} and \code{dfbetas()}, a matrix with one row per cluster
#' or observation, in the same order, and one column per coefficient.
#'
#' @references
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' Preisser, J.S. and Qaqish, B.F. (1996) Deletion diagnostics for generalised
#' estimating equations. \emph{Biometrika}, \bold{83}, 551--562.
#'
#' Vanegas, L.H., Rondon, L.M. and Paula, G.A. (2023) Generalized Estimating
#' Equations using the new R package glmtoolbox. \emph{The R Journal},
#' \bold{15}, 105--133.
#'
#' @seealso \code{\link{local_influence}}, \code{\link{residuals.geer}},
#'   \code{\link{vcov.geer}}, \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   data = epilepsy,
#'   id = id,
#'   repeated = visit,
#'   family = poisson(link = "log"),
#'   corstr = "exchangeable"
#' )
#' head(hatvalues(fit))
#' head(dfbeta(fit))
#' head(dfbetas(fit, cov_type = "robust"))
#' head(cooks.distance(fit))
#' head(cooks.distance(fit, level = "observations"))
#' head(cooks.distance(fit, parm = "treatmentprogabide"))
#'
#' @export
hatvalues.geer <- function(model,
                           level = c("clusters", "observations"),
                           ...) {
  check_unused_dots(list(...), "hatvalues.geer")
  model <- check_geer_object(model, name = "model")
  level <- match_influence_level(level[1L])
  compute_influence_hatvalues(model, level)
}


#' @rdname influence.geer
#' @method dfbeta geer
#' @export
dfbeta.geer <- function(model,
                        level = c("clusters", "observations"),
                        method = c("one-step", "exact"),
                        ...) {
  check_unused_dots(list(...), "dfbeta.geer")
  model <- check_geer_object(model, name = "model")
  level <- match_influence_level(level[1L])
  method <- match_influence_method(method[1L])
  compute_influence_dfbeta(model, level, method, "dfbeta.geer")
}


#' @rdname influence.geer
#' @method dfbetas geer
#' @export
dfbetas.geer <- function(model,
                         level = c("clusters", "observations"),
                         method = c("one-step", "exact"),
                         cov_type = "naive",
                         ...) {
  check_unused_dots(list(...), "dfbetas.geer")
  model <- check_geer_object(model, name = "model")
  level <- match_influence_level(level[1L])
  method <- match_influence_method(method[1L])
  cov_type <- match_cov_type(cov_type, default = "naive")
  compute_influence_dfbetas(model, level, method, cov_type, "dfbetas.geer")
}


#' @rdname influence.geer
#' @method cooks.distance geer
#' @export
cooks.distance.geer <- function(model,
                                level = c("clusters", "observations"),
                                method = c("one-step", "exact"),
                                cov_type = "naive",
                                parm,
                                ...) {
  check_unused_dots(list(...), "cooks.distance.geer")
  model <- check_geer_object(model, name = "model")
  level <- match_influence_level(level[1L])
  method <- match_influence_method(method[1L])
  cov_type <- match_cov_type(cov_type, default = "naive")
  parm <- match_parm(
    if (missing(parm)) NULL else parm,
    names(model$coefficients)
  )
  compute_cooks_distance(
    model,
    level = level,
    method = method,
    cov_type = cov_type,
    parm = parm,
    context = "cooks.distance.geer"
  )
}
