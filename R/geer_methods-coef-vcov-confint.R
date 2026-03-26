#' @title
#' Extract Model Coefficients from a \code{geer} Object
#'
#' @aliases coef.geer coef coefficients
#' @method coef geer
#'
#' @description
#' Extracts model coefficients from a \code{geer} object. The function
#' \code{coefficients} is an alias.
#'
#' @inheritParams add1.geer
#'
#' @return
#' A named numeric vector of estimated regression coefficients.
#'
#' @examples
#' data("leprosy", package = "geer")
#' fit <- geewa(
#'   formula = bacilli ~ factor(period) + factor(period):treatment,
#'   family = poisson(link = "log"),
#'   id = id,
#'   data = leprosy
#' )
#' coef(fit)
#' names(coef(fit))
#'
#' data("cerebrovascular", package = "geer")
#' fit2 <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   id = id,
#'   data = cerebrovascular
#' )
#' coef(fit2)
#'
#' @export
coef.geer <- function(object, ...) {
  object <- check_geer_object(object)
  object$coefficients
}



#' @title
#' Confidence Intervals for Model Parameters from a \code{geer} Object
#'
#' @aliases confint confint.geer
#' @method confint geer
#'
#' @description
#' Compute Wald-type confidence intervals for one or more regression parameters
#' from a fitted \code{geer} model.
#'
#' @inheritParams add1.geer
#' @inheritParams stats::confint
#'
#' @details
#' Confidence intervals are computed as
#' \eqn{\hat\beta \pm z_{1-\alpha/2}\,\mathrm{SE}(\hat\beta)},
#' where standard errors are obtained from \code{vcov(object, cov_type = cov_type)}.
#' The resulting intervals rely on the usual large-sample normal approximation.
#'
#' @inherit stats::confint.default return
#'
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   id = id,
#'   data = cerebrovascular,
#'   orstr = "exchangeable"
#' )
#' confint(fit)
#' confint(fit, parm = "treatmentactive")
#' confint(fit, cov_type = "naive")
#' @export
confint.geer <- function(object,
                         parm,
                         level = 0.95,
                         cov_type = c("robust", "bias-corrected",
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
  alpha <- (1 - level) / 2
  alpha <- c(alpha, 1 - alpha)
  pct <- format_perc(alpha, 3)
  percentiles <- qnorm(alpha)
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

#' @title
#' Extract Variance-Covariance Matrix from a \code{geer} Object
#'
#' @aliases vcov vcov.geer
#' @method vcov geer
#'
#' @description
#' Extract the variance–covariance matrix of the regression parameters from a
#' fitted \code{geer} object. The parameters correspond to those returned by
#' \code{\link{coef}}.
#'
#' @inheritParams add1.geer
#'
#' @details
#' The form of the covariance estimator is controlled by the argument
#' \code{cov_type}:
#' \itemize{
#'   \item \code{"robust"} -- the sandwich (robust) covariance estimator
#'     \cite{Liang and Zeger (1986)}
#'   \item \code{"naive"} -- the model-based covariance estimator
#'     \cite{Liang and Zeger (1986)}
#'   \item \code{"df-adjusted"} -- the small-sample adjusted covariance matrix
#'     \cite{MacKinnon (1985)}
#'   \item \code{"bias-corrected"} -- the bias-corrected covariance estimator
#'     \cite{Morel et al. (2003)}
#' }
#'
#' @return
#' A square numeric matrix of estimated covariances between regression
#' coefficients. Rows and columns are named according to the coefficient names
#' returned by \code{\link[stats]{coef}}.
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using generalized
#' linear models \emph{Biometrika} \bold{73}, 13-–22.
#'
#' MacKinnon, J.G. (1985) Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics} \bold{29}, 305–-325.
#'
#' Morel, J.G., Bokossa, M.C. and Neerchal N.K. (2003) Small sample correction
#' for the variance of GEE estimators. \emph{Biometrical Journal} \bold{45},
#' 395–-409.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fitted_model <- geewa_binary(formula = ecg ~ period + treatment,
#'                              id = id,
#'                              data = cerebrovascular,
#'                              link = "logit",
#'                              orstr = "exchangeable")
#' vcov(fitted_model)
#' vcov(fitted_model, cov_type = "bias-corrected")
#'
#' @export
vcov.geer <- function(object,
                      cov_type = c("robust", "bias-corrected",
                                   "df-adjusted", "naive"),
                      ...) {
  object <- check_geer_object(object)
  cov_type <- match.arg(cov_type)
  switch(
    cov_type,
    robust = object$robust_covariance,
    naive = object$naive_covariance,
    `bias-corrected` = object$bias_corrected_covariance,
    `df-adjusted` = {
      sample_size <- object$clusters_no
      p <- ncol(object$robust_covariance)
      if (sample_size <= p) {
        stop("clusters_no must be > number of coefficients", call. = FALSE)
      }
      (sample_size / (sample_size - p)) * object$robust_covariance
    }
  )
}
