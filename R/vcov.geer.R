#' @title
#' Calculate variance-covariance matrix from \code{geer} objects
#'
#' @aliases vcov vcov.geer
#'
#' @method vcov geer
#'
#' @description
#' Returns the variance-covariance matrix of the main parameters of a fitted
#' model geer object. The "main" parameters of model correspond to those returned
#' by [stats::coef()].
#'
#' @param object a fitted model \code{geer} object.
#' @param cov_type character indicating whether the sandwich (robust) covariance
#'        matrix (\code{cov_type = "robust"}), the model-based (naive) covariance
#'        matrix (\code{cov_type = "naive"}), the bias-corrected covariance
#'        matrix (\code{cov_type = "bias-corrected"}) or the degrees of freedom adjusted
#'        covariance matrix (\code{cov_type = "df-adjusted"}) should be returned. By
#'        default, the robust covariance matrix is returned.
#' @param ... additional argument(s) for methods.
#'
#' @details
#' If \code{cov_type = "robust"}, then the so-called sandwich or
#' robust covariance matrix is returned, see \cite{Liang and Zeger (1986)}.
#'
#' If \code{cov_type = "naive"}, then the so-called naive or
#' model-based covariance matrix is returned, see \cite{Liang and Zeger (1986)}.
#'
#' If \code{cov_type = "df-adjusted"}, then the adjusted covariance matrix proposed
#' by \cite{MacKinnon (1985)} is returned.
#'
#' If \code{cov_type = "bias-corrected"}, then the bias-corrected covariance matrix
#' proposed by \cite{Morel, Bokossa and Neerchal (2003)} is returned.
#'
#' @return
#' A matrix of the estimated covariances between the parameter estimates
#' in the linear predictor of the model. This should have row and column
#' names corresponding to the parameter names given by the [stats::coef()] method.
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) A comparison of two bias-corrected covariance
#' estimators for generalized estimating equations. \emph{Biometrika} \bold{73},
#' 13-–22.
#'
#' MacKinnon, J.G. (1985) Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics} \bold{29}, 305–-325.
#'
#' Morel J.G., Bokossa M.C. and Neerchal N.K, (2003) Small sample correction
#' for the variance of GEE estimators. \emph{Biometrical Journal} \bold{45},
#' 395–-409.
#'
#' @examples
#' data("cerebrovascular")
#' fitted_model <- geewa_binary(formula = ecg ~ period + treatment,
#'                              id = id,
#'                              data = cerebrovascular,
#'                              link = "logit",
#'                              orstr = "exchangeable")
#' vcov(fitted_model)
#' vcov(fitted_model, cov_type = "bias-corrected")
#'
#' @export

vcov.geer <- function(object, cov_type = "robust", ...) {
  icheck <- pmatch(cov_type,
                   c("robust", "naive", "df-adjusted", "bias-corrected"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown method for the covariance matrix")
  if (cov_type == "robust") {
    ans <- object$robust_covariance
    } else if (cov_type == "bias-corrected") {
      total_obs <- object$obs_no
      sample_size <- object$clusters_no
      parameters_no <- length(object$coefficients)
      kappa <-
        ((total_obs - 1)/(total_obs - parameters_no)) *
        (sample_size / (sample_size - 1))
      delta <-
        min(parameters_no / (sample_size - parameters_no),
            0.5)
      robust_covariance <- object$robust_covariance
      naive_covariance <- object$naive_covariance
      trace_robust_naive_inverse <- sum(diag(solve(naive_covariance, robust_covariance)))
      ksi <-
        max(1,
            trace_robust_naive_inverse/parameters_no)
      ans <- kappa * robust_covariance + delta * ksi * naive_covariance
    } else if (cov_type == "df-adjusted") {
      sample_size <- object$clusters_no
      parameters_no <- length(object$coefficients)
      ans <-
        (sample_size / (sample_size - parameters_no)) *
        object$robust_covariance
    } else {
      ans <- object$naive_covariance
    }
  ans
}
