#' @title
#' Extract Variance-Covariance Matrix from a geer Object
#'
#' @aliases vcov vcov.geer
#' @method vcov geer
#'
#' @description
#' Extracts the variance-covariance matrix of the estimated regression
#' parameters from a fitted \code{geer} object.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for inference on the regression parameters. Options are the bias-corrected estimator
#'   (\code{"bias-corrected"}), the sandwich or robust estimator (\code{"robust"}), the degrees-of-freedom adjusted estimator
#'   (\code{"df-adjusted"}), and the model-based or naive estimator
#'   (\code{"naive"}). Defaults to \code{"bias-corrected"}.
#' @param ... additional arguments passed to or from other methods.
#'
#' @details
#' The form of the covariance estimator is controlled by \code{cov_type}:
#' \describe{
#'   \item{\code{"bias-corrected"}}{the bias-corrected covariance estimator
#'   (Morel et al., 2003).}
#'   \item{\code{"robust"}}{the sandwich (robust) covariance estimator
#'   (Liang and Zeger, 1986).}
#'   \item{\code{"df-adjusted"}}{the degrees-of-freedom adjusted covariance
#'   estimator (MacKinnon, 1985).}
#'   \item{\code{"naive"}}{the model-based covariance estimator
#'   (Liang and Zeger, 1986).}
#' }
#'
#' @return
#' A square numeric matrix of estimated covariances between regression
#' coefficients. Rows and columns are named according to the coefficient names
#' returned by \code{\link{coef.geer}}.
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, \bold{73}, 13--22.
#'
#' MacKinnon, J.G. (1985) Some heteroskedasticity-consistent covariance matrix
#' estimators with improved finite sample properties. \emph{Journal of
#' Econometrics}, \bold{29}, 305--325.
#'
#' Morel, J.G., Bokossa, M.C. and Neerchal, N.K. (2003) Small sample
#' correction for the variance of GEE estimators. \emph{Biometrical Journal},
#' \bold{45}, 395--409.
#'
#' @seealso \code{\link{coef.geer}}, \code{\link{confint.geer}},
#'   \code{\link{summary.geer}}, \code{\link{tidy.geer}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' vcov(fit)
#' vcov(fit, cov_type = "robust")
#'
#' @export
vcov.geer <- function(object,
                      cov_type = c("bias-corrected", "robust",
                                   "df-adjusted", "naive"),
                      ...) {
  object <- check_geer_object(object)
  cov_type <- match.arg(cov_type)
  switch(
    cov_type,
    robust = object$robust_covariance,
    naive = object$naive_covariance,
    `bias-corrected` = object$bias_corrected_covariance,
    `df-adjusted` = compute_df_adjusted_covariance(
      robust_covariance = object$robust_covariance,
      clusters_no = object$clusters_no,
      coef_no = ncol(object$robust_covariance),
      context = "vcov"
    )
  )
}


#' @title
#' Extract Model Coefficients from a geer Object
#'
#' @aliases coef.geer coef coefficients
#' @method coef geer
#'
#' @description
#' Extracts the estimated regression coefficients from a fitted \code{geer}
#' object. \code{coefficients} is an alias for \code{coef}.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param ... additional arguments passed to or from other methods.
#'
#' @return
#' A named numeric vector of estimated regression coefficients. The names
#' correspond to the columns of the model matrix.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{confint.geer}},
#'   \code{\link{summary.geer}}.
#'
#' @examples
#' data("leprosy", package = "geer")
#' fit <- geewa(
#'   formula = bacilli ~ factor(period) + factor(period):treatment,
#'   family = poisson(link = "log"),
#'   data = leprosy,
#'   id = id
#' )
#' coef(fit)
#'
#' data("cerebrovascular", package = "geer")
#' fit_bin <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id
#' )
#' coef(fit_bin)
#'
#' @export
coef.geer <- function(object, ...) {
  object <- check_geer_object(object)
  object$coefficients
}


#' @title
#' Confidence Intervals for Model Parameters from a geer Object
#'
#' @aliases confint confint.geer
#' @method confint geer
#'
#' @description
#' Computes Wald-type confidence intervals for one or more regression
#' parameters from a fitted \code{geer} object.
#'
#' @inheritParams vcov.geer
#' @inheritParams stats::confint
#'
#' @details
#' Confidence intervals are computed as
#' \eqn{\hat{\beta} \pm z_{1-\alpha/2} \, \mathrm{SE}(\hat{\beta})},
#' where standard errors are obtained from
#' \code{vcov(object, cov_type = cov_type)}. The covariance estimator used is
#' controlled by \code{cov_type}; see \code{\link{vcov.geer}} for details.
#' The resulting intervals rely on the usual large-sample normal approximation.
#'
#' @return
#' A matrix with columns giving lower and upper confidence limits for each
#' parameter. The columns are labeled by the corresponding tail probabilities
#' in percent, for example \code{"2.5\%"} and \code{"97.5\%"} when
#' \code{level = 0.95}.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{coef.geer}},
#'   \code{\link{tidy.geer}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' confint(fit)
#' confint(fit, parm = "treatmentactive")
#' confint(fit, cov_type = "naive")
#'
#' @export
confint.geer <- function(object,
                         parm,
                         level = 0.95,
                         cov_type = c("bias-corrected", "robust",
                                      "df-adjusted", "naive"),
                         ...) {
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be a single number in (0, 1)", call. = FALSE)
  }
  object <- check_geer_object(object)
  cov_type <- match.arg(cov_type)
  beta <- coef(object)
  beta_names <- names(beta)
  if (missing(parm)) {
    parm <- beta_names
  } else if (is.numeric(parm)) {
    parm <- beta_names[parm]
  } else if (!is.character(parm)) {
    stop("'parm' must be missing, numeric, or character", call. = FALSE)
  }
  if (anyNA(parm) || !all(parm %in% beta_names)) {
    stop("invalid 'parm' specification", call. = FALSE)
  }
  conf_probs <- (1 - level) / 2
  conf_probs <- c(conf_probs, 1 - conf_probs)
  pct <- format_percent(conf_probs, 3)
  percentiles <- qnorm(conf_probs)
  ans <- matrix(
    NA_real_,
    nrow = length(parm),
    ncol = 2L,
    dimnames = list(parm, pct)
  )
  vcov_matrix <- vcov(object, cov_type = cov_type)
  standard_errors <- sqrt(pmax(diag(vcov_matrix), 0))[parm]
  ans[] <- beta[parm] + standard_errors %o% percentiles
  ans
}
