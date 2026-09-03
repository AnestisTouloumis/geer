#' @title
#' Bias-Reduced and Penalized Generalized Estimating Equations
#'
#' @description
#' Fits marginal models for repeated or clustered responses using
#' Generalized Estimating Equations (GEE). Supported estimation methods include
#' the traditional GEE, bias-reducing GEE, bias-corrected GEE, and
#' Jeffreys-type penalized GEE. Continuous, binary and count responses are handled by
#' \code{\link{geewa}}, while binary responses can also be handled by
#' \code{\link{geewa_binary}} through an odds-ratio parameterization.
#'
#' @author Anestis Touloumis \email{A.Touloumis@@brighton.ac.uk}
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, \bold{73}, 13--22.
#'
#' Touloumis, A. (2026) Jeffreys-Type Penalized GEE for Correlated
#' Binary Data with an Odds-Ratio Parameterization. \emph{Preprint}.
#' \url{https://arxiv.org/abs/2606.16058}
#'
#' Touloumis, A. (2026) Bias-Reduced GEE via Adjusted Estimating Equations, with Odds-Ratio Extensions.
#' \emph{Preprint}. \url{https://arxiv.org/abs/2606.16043}
#'
#' @seealso
#' Main functions:
#' \itemize{
#'   \item \code{\link{geewa}} for continuous, binary and count responses.
#'   \item \code{\link{geewa_binary}} for binary responses via an odds-ratio
#'   parameterization (preferred).
#'   \item \code{\link{geer_control}} for convergence and fitting options.
#'   \item \code{\link{geecriteria}} for model selection criteria.
#'   \item \code{\link{summary.geer}}, \code{\link{tidy.geer}}, and
#'   \code{\link{glance.geer}} for model summaries.
#'   \item \code{\link{anova.geer}}, \code{\link{add1.geer}},
#'   \code{\link{drop1.geer}}, and \code{\link{step_p}} for model building.
#'   \item \code{\link{vcov.geer}}, \code{\link{confint.geer}},
#'   \code{\link{predict.geer}}, \code{\link{fitted.geer}},
#'   \code{\link{residuals.geer}}, \code{\link{runs_test}},
#'   \code{\link{mcar_little_test}}, \code{\link{mcar_homoscedasticity_test}},
#'   \code{\link{mcar_logistic_test}}, and
#'   \code{\link{frechet_bounds_cor}}
#'   for inference and model diagnostics.
#'   \item \code{\link{marginaleffects-support}} for post-estimation support
#'   through the \pkg{marginaleffects} package.
#' }
#'
#' @name geer-package
#'
#' @useDynLib geer, .registration = TRUE
#' @import Rcpp
#' @importFrom brglm2 brglmFit brglm_control
#' @importFrom generics tidy glance
## Two groups are imported rather than called with an explicit stats::
## prefix. First, the generics for which geer registers S3 methods, which
## must be resolvable when the S3method() directives are processed. Second,
## the family constructors, which are reached indirectly by name through
## match.fun() in normalize_family() and do.call() in geewa(), so a call-site
## prefix cannot be used. Every other stats/utils function is called with an
## explicit stats:: or utils:: prefix.
#' @importFrom stats add1 anova coef confint drop1 fitted model.matrix
#' @importFrom stats predict residuals vcov
#' @importFrom stats Gamma binomial gaussian inverse.gaussian poisson
#' @importFrom stats quasi quasibinomial quasipoisson
#' @keywords package
"_PACKAGE"
