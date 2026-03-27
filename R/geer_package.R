#' GEE Solver in R
#'
#' Solving generalized estimating equations while implementing bias-reducing,
#' bias-corrected and penalized estimators.
#'
#' @name geer-package
#' @aliases geer
#' @author Anestis Touloumis
#'
#' Maintainer: Anestis Touloumis \email{A.Touloumis@@brighton.ac.uk}
#'
#' @useDynLib geer, .registration = TRUE
#' @import Rcpp
#' @importFrom brglm2 brglmFit brglm_control
#' @importFrom methods missingArg
#' @importFrom stats ave .getXlevels add.scope add1 as.formula binomial coef
#' @importFrom stats coefficients dbinom delete.response dnorm dpois drop.scope
#' @importFrom stats drop1 factor.scope family formula gaussian glm glm.fit
#' @importFrom stats model.extract model.frame model.matrix model.offset
#' @importFrom stats model.response model.weights pchisq pnorm printCoefmat
#' @importFrom stats qlogis qnorm terms update update.formula vcov
#' @importFrom stats Gamma inverse.gaussian poisson confint na.pass
#' @importFrom utils combn
#' @importFrom generics tidy
#' @export
generics::tidy
#' @importFrom generics glance
#' @export
generics::glance
#' @keywords internal
#'@importFrom generics tidy
"_PACKAGE"
