#' @title
#' Bias-Reduced and Penalized Generalized Estimating Equations
#'
#' @description
#' Fits marginal models for clustered or repeated responses by generalized
#' estimating equations (GEE; Liang and Zeger, 1986), using ordinary,
#' bias-reduced, or bias-corrected GEE (Touloumis, 2026a) or Jeffreys-type
#' penalized GEE (Touloumis, 2026b). \code{\link{geewa}} handles continuous,
#' binary, and count responses with a working correlation structure, and
#' \code{\link{geewa_binary}} handles binary responses with a working
#' odds-ratio structure.
#'
#' @author Anestis Touloumis \email{A.Touloumis@@brighton.ac.uk}
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, \bold{73}, 13--22.
#'
#' Touloumis, A. (2026a) Bias-reduced GEE via adjusted estimating equations,
#' with odds-ratio extensions. \emph{arXiv preprint}.
#' \doi{10.48550/arXiv.2606.16043}
#'
#' Touloumis, A. (2026b) Jeffreys-type penalized GEE for correlated binary data
#' with an odds-ratio parameterization. \emph{arXiv preprint}.
#' \doi{10.48550/arXiv.2606.16058}
#'
#' @seealso
#' \itemize{
#'   \item Model fitting: \code{\link{geewa}}, \code{\link{geewa_binary}},
#'     \code{\link{geer_control}}.
#'   \item Summaries and inference: \code{\link{summary.geer}},
#'     \code{\link{vcov.geer}}, \code{\link{confint.geer}},
#'     \code{\link{predict.geer}}, \code{\link{tidy.geer}},
#'     \code{\link{glance.geer}}.
#'   \item Model building: \code{\link{anova.geer}}, \code{\link{add1.geer}},
#'     \code{\link{step_p}}, \code{\link{step_criterion}},
#'     \code{\link{geecriteria}}.
#'   \item Diagnostics: \code{\link{residuals.geer}}, \code{\link{runs_test}},
#'     \code{\link{frechet_bounds_cor}}, \code{\link{mcar_little_test}},
#'     \code{\link{mcar_logistic_test}},
#'     \code{\link{mcar_homoscedasticity_test}}.
#'   \item Other packages: \code{\link{emmeans-support}},
#'     \code{\link{marginaleffects-support}}.
#' }
#'
#' @name geer-package
#'
#' @useDynLib geer, .registration = TRUE
#' @import Rcpp
#' @importFrom brglm2 brglmFit brglm_control
#' @importFrom generics tidy glance
#' @importFrom stats add1 anova coef confint drop1 fitted model.matrix
#' @importFrom stats predict residuals vcov
#' @importFrom stats Gamma binomial gaussian inverse.gaussian poisson
#' @importFrom stats quasi quasibinomial quasipoisson
#' @keywords package
"_PACKAGE"
