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
#'
#' @importFrom stats gaussian glm model.extract model.matrix model.response
#' @importFrom stats coef pnorm printCoefmat binomial qnorm vcov coefficients
#' @importFrom stats update pchisq formula family .getXlevels delete.response
#' @importFrom stats model.frame model.offset as.formula update.formula
#' @importFrom stats drop1 add1 add.scope drop.scope terms formula factor.scope
#' @importFrom stats glm.fit model.weights
#' @importFrom methods missingArg
#' @importFrom utils combn
#' @keywords internal
"_PACKAGE"
