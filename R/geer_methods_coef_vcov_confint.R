#' @title
#' Extract Variance-Covariance Matrix from a geer Object
#'
#' @method vcov geer
#'
#' @description
#' Extracts the variance-covariance matrix of the estimated regression
#' parameters from a fitted \code{geer} object.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for inference: \code{"bias-corrected"}, \code{"robust"},
#'   \code{"df-adjusted"}, \code{"jackknife"}, or \code{"naive"}. See Details.
#'   Defaults to \code{"bias-corrected"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @details
#' The covariance estimator is selected by \code{cov_type}:
#' \describe{
#'   \item{\code{"bias-corrected"}}{the bias-corrected sandwich estimator of
#'     Morel, Bokossa and Neerchal (2003).}
#'   \item{\code{"robust"}}{the sandwich estimator of Liang and Zeger (1986).}
#'   \item{\code{"df-adjusted"}}{the sandwich estimator with the
#'     degrees-of-freedom adjustment of MacKinnon and White (1985).}
#'   \item{\code{"jackknife"}}{the leave-one-cluster-out jackknife estimator.}
#'   \item{\code{"naive"}}{the model-based estimator of Liang and Zeger
#'     (1986).}
#' }
#' Both the bias-corrected and the degrees-of-freedom adjusted estimators
#' require more clusters than regression parameters: otherwise the former is
#' \code{NA} and the latter gives an error.
#'
#' With \eqn{K} clusters and \eqn{\hat\beta_{(i)}}{beta_(i)} the estimate obtained after
#' deleting cluster \eqn{i}, the jackknife estimator is
#' \deqn{\frac{K - 1}{K} \sum_{i=1}^K (\hat\beta_{(i)} - \bar\beta)
#' (\hat\beta_{(i)} - \bar\beta)^\top, \qquad
#' \bar\beta = \frac{1}{K}\sum_{i=1}^K \hat\beta_{(i)}.}{(K - 1)/K sum_i (beta_(i) - betabar)(beta_(i) - betabar)',   betabar = (1/K) sum_i beta_(i).}
#' Each \eqn{\hat\beta_{(i)}}{beta_(i)} is a full refit by the original estimation
#' method, not a one-step approximation, with the working association
#' parameters held at their full-data estimates. The dispersion is
#' re-estimated if it was estimated in the original fit. For the two-stage
#' methods (\code{"bcgee-*"}, \code{"opgee-jeffreys"}, and
#' \code{"hpgee-jeffreys"}) the first stage is repeated in each refit. At
#' least two clusters are required, and an error is signaled if a refit does
#' not converge. The jackknife estimate is cached in the fitted object.
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
#' MacKinnon, J.G. and White, H. (1985) Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics}, \bold{29}, 305--325.
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
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' vcov(fit)
#' vcov(fit, cov_type = "robust")
#'
#' @export
vcov.geer <- function(object,
                      cov_type = c("bias-corrected", "robust", "df-adjusted",
                                   "jackknife", "naive"),
                      ...) {
  check_unused_dots(list(...), "vcov.geer")
  object <- check_geer_object(object)
  cov_type <- match.arg(cov_type)
  switch(
    cov_type,
    robust = object$robust_covariance,
    naive = object$naive_covariance,
    `bias-corrected` = object$bias_corrected_covariance,
    jackknife = get_cached_jackknife_covariance(object),
    `df-adjusted` = compute_df_adjusted_covariance(
      robust_covariance = object$robust_covariance,
      clusters_no = object$clusters_no,
      coef_no = object$rank,
      context = "vcov"
    )
  )
}


#' @title
#' Extract Model Coefficients from a geer Object
#'
#' @method coef geer
#'
#' @description
#' Extracts the estimated regression coefficients from a fitted \code{geer}
#' object. \code{coefficients} is an alias for \code{coef}.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @return
#' A named numeric vector of estimated regression coefficients. The names
#' correspond to the columns of the model matrix.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{confint.geer}},
#'   \code{\link{summary.geer}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' coef(fit)
#'
#' @export
coef.geer <- function(object, ...) {
  check_unused_dots(list(...), "coef.geer")
  object <- check_geer_object(object)
  object$coefficients
}


#' @title
#' Confidence Intervals for Model Parameters from a geer Object
#'
#' @method confint geer
#'
#' @description
#' Computes Wald-type confidence intervals for one or more regression
#' parameters from a fitted \code{geer} object.
#'
#' @inheritParams vcov.geer
#' @inheritParams stats::confint
#' @param level a single number strictly between 0 and 1 specifying the
#'   confidence level. Defaults to \code{0.95}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for the standard errors: \code{"bias-corrected"}, \code{"robust"},
#'   \code{"df-adjusted"}, \code{"jackknife"}, or \code{"naive"}. See
#'   \code{\link{vcov.geer}}. Defaults to \code{"bias-corrected"}.
#'
#' @details
#' The intervals are
#' \eqn{\hat\beta \pm z_{1-\alpha/2}\,\mathrm{SE}(\hat\beta)}{beta +/- z_(1 - alpha/2) SE(beta)}, with standard
#' errors from \code{vcov(object, cov_type = cov_type)}.
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
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' confint(fit)
#' confint(fit, parm = "treatmentactive", cov_type = "robust")
#'
#' @export
confint.geer <- function(object,
                         parm,
                         level = 0.95,
                         cov_type = c("bias-corrected", "robust",
                                      "df-adjusted", "jackknife", "naive"),
                         ...) {
  check_unused_dots(list(...), "confint.geer")
  object <- check_geer_object(object)
  check_probability_open(level, "level")
  cov_type <- match.arg(cov_type)
  beta <- object$coefficients
  beta_names <- names(beta)
  if (missing(parm)) {
    parm <- beta_names
  } else if (is.numeric(parm)) {
    parm <- beta_names[parm]
  } else if (!is.character(parm)) {
    stop("'parm' must be missing, numeric, or character", call. = FALSE)
  }
  if (anyNA(parm) || !all(parm %in% beta_names)) {
    stop("'parm' must index or name model coefficients", call. = FALSE)
  }
  conf_probs <- (1 - level) / 2
  conf_probs <- c(conf_probs, 1 - conf_probs)
  pct <- format_percent(conf_probs, 3)
  percentiles <- stats::qnorm(conf_probs)
  ans <- matrix(
    NA_real_,
    nrow = length(parm),
    ncol = 2L,
    dimnames = list(parm, pct)
  )
  vcov_matrix <- stats::vcov(object, cov_type = cov_type)
  standard_errors <- compute_standard_errors_or_stop(
    vcov_matrix,
    parm = parm,
    context = "confint.geer"
  )
  ans[] <- beta[parm] + standard_errors %o% percentiles
  ans
}
