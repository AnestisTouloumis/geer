#' GEE Solver in R
#'
#' Solving generalized estimating equation while implementing bias-reducing,
#' bias-corrected and penalized estimators.
#'
#' @name geer-package
#' @aliases geer
#' @docType package
#' @author Anestis Touloumis
#'
#' Maintainer: Anestis Touloumis \email{A.Touloumis@@brighton.ac.uk}
#'
#' @useDynLib geer, .registration = TRUE
#' @import Rcpp
#' @importFrom brglm2 brglmFit
#'
#' @importFrom stats gaussian glm model.extract model.matrix model.response
#' @importFrom stats coef pnorm printCoefmat binomial qnorm vcov coefficients
#' @importFrom stats update pchisq formula family
#' @importFrom utils combn
#' @keywords package
NULL
