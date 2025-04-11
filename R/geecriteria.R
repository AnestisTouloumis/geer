#' @title
#' Variable and Covariance Selection Criteria for \code{geer} Objects
#'
#' @description
#' Reports commonly used criteria for variable selection and association structure
#' selection for one or several fitted models.
#'
#' @return
#' A data frame with the QIC, CIC, RJC, QICu, GESSC, GPC and the number
#' of regression parameters (including intercepts).
#'
#' @inheritParams anova.geer
#' @param cov_type character indicating the type of covariance matrix estimator
#'        which should be used to the calculation of the criteria. Options
#'        include the sandwich or robust estimator (\code{"robust"}), the
#'        bias-corrected estimator (\code{"bias-corrected"}), the degrees of
#'        freedom adjusted estimator (\code{"df-adjusted"}) and the model-based
#'        or naive estimator (\code{"naive"}). By default, the robust covariance
#'        estimator is used.
#' @param digits integer indicating the number of decimal places for the GEE
#'        criteria to be rounded at. By default, \code{digits = 2}.
#'
#' @details
#' The Quasi Information Criterion (QIC), the Correlation Information Criterion
#' (CIC), the Rotnitzky and Jewell Criterion (RJC), the Generalized Error Sum of
#' Squares (GESSC) and the Gaussian Pseudolikelihood Criterion (GPC) are used for
#' selecting the best association structure.
#' The QICu criterion is used for selecting the best subset of covariates.
#'
#' When choosing between ordinary, bias-reducing and penalized GEE models with
#' the same subset of covariates, the model with the smallest value of QIC, CIC,
#' RJC or GESSC or the model with the higher GPC should be preferred. When
#' choosing between ordinary, bias-reducing and penalized models with different
#' number of covariates, the model with the smallest QICu value should be
#' preferred.
#'
#' @references
#' Carey, V.J. and Wang, Y.G. (2011) Working covariance model selection for
#' generalized estimating equations. \emph{Statistics in Medicine}, \bold{30},
#' 3117--3124.
#'
#' Chaganty, N.R. and Shults, J. (1999) On eliminating the asymptotic bias in
#' the quasi-least squares estimate of the correlation parameter.
#' \emph{Journal of Statistical Planning and Inference}, \bold{76}, 145--161.
#'
#' Hin, L.Y. and Wang, Y.G. (2009) Working correlation structure
#' identification in generalized estimating equations. \emph{Statistics in
#' Medicine} \bold{28}, 642--658.
#'
#' Pan, W. (2001) Akaike's information criterion in generalized
#' estimating equations. \emph{Biometrics} \bold{57}, 120--125.
#'
#' Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
#' parameters in semiparametric generalized linear models for cluster correlated
#' data. \emph{Biometrika} \bold{77}, 485--497.
#'
#' @examples
#' data("cerebrovascular")
#' fitted_gee <- geewa_binary(formula = ecg ~ period * treatment,
#'                            id = id,
#'                            data = cerebrovascular,
#'                            link = "logit",
#'                            orstr = "exchangeable",
#'                            method = "gee")
#' fitted_brgee <- update(fitted_gee, method = "brgee_robust")
#' geecriteria(fitted_gee, fitted_brgee)
#'
#' @export
geecriteria <- function(object, ..., cov_type = "robust", digits = 2) {
  if (length(list(...))) {
    all_models <- list(object, ...)
    if (any(lapply(all_models, function(x) class(x)[1]) != "geer"))
      stop("Only 'geer' objects are supported")
    results <- lapply(all_models, compute_criteria, cov_type, digits)
    check <- sapply(all_models, function(x) length(x$obs_no))
    if (any(check != check[1]))
      warning("models do not have the same number of observations")
    ans <- do.call("rbind", results)
    Call <- match.call()
    remove_names <- c(which(names(Call) == "cov_type"),
                      which(names(Call) == "digits"))
    row.names(ans) <- as.character(Call[-c(1, remove_names)])
  } else {
    ans <- compute_criteria(object, cov_type, digits)
  }
  ans
}
