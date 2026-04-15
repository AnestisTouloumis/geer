#' Generalized Estimating Equations in R
#'
#' @description
#' Fits marginal models for independent, repeated, or clustered responses using
#' Generalized Estimating Equations (GEE). Supported estimation methods include
#' the traditional GEE, bias-reducing GEE, bias-corrected GEE, and
#' Jeffreys-prior penalized GEE. Continuous and count responses are handled by
#' \code{\link{geewa}}, while binary responses are handled by
#' \code{\link{geewa_binary}} through an odds-ratio parameterization.
#'
#' @author Anestis Touloumis \email{A.Touloumis@@brighton.ac.uk}
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, \bold{73}, 13--22.
#'
#' Touloumis, A. (2026) Bias reduction in generalized estimating equations.
#' Preprint.
#'
#' Touloumis, A. (2026) Jeffreys-prior penalized GEE for correlated binary
#' data with an odds-ratio parameterization. Preprint.
#'
#' @seealso
#' Main functions:
#' \itemize{
#'   \item \code{\link{geewa}} for continuous and count responses.
#'   \item \code{\link{geewa_binary}} for binary responses via an odds-ratio
#'   parameterization.
#'   \item \code{\link{geer_control}} for convergence and fitting options.
#'   \item \code{\link{geecriteria}} for model selection criteria.
#'   \item \code{\link{summary.geer}}, \code{\link{tidy.geer}}, and
#'   \code{\link{glance.geer}} for model summaries.
#'   \item \code{\link{anova.geer}}, \code{\link{add1.geer}},
#'   \code{\link{drop1.geer}}, and \code{\link{step_p}} for model building.
#' }
#'
#' @name geer-package
#'
#' @useDynLib geer, .registration = TRUE
#' @import Rcpp
#' @importFrom brglm2 brglmFit brglm_control
#' @importFrom stats ave .getXlevels add.scope add1 as.formula binomial coef
#' @importFrom stats coefficients dbinom delete.response dnorm dpois drop.scope
#' @importFrom stats drop1 factor.scope family formula gaussian glm glm.fit
#' @importFrom stats model.extract model.frame model.matrix model.offset
#' @importFrom stats model.response model.weights pchisq pnorm printCoefmat
#' @importFrom stats qlogis qnorm terms update update.formula vcov
#' @importFrom stats Gamma inverse.gaussian poisson confint na.pass
#' @importFrom utils combn
#' @importFrom generics tidy glance
#' @keywords package
"_PACKAGE"
