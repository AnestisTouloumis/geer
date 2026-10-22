#' @title
#' Little's Test of Missing Completely at Random
#'
#' @description
#' Performs Little's (1988) test of the hypothesis that the missing values are
#' missing completely at random (MCAR).
#'
#' @param object a fitted model object of class \code{"geer"}, or a numeric
#'   matrix or data frame with one row per independent unit and one column per
#'   variable or repeated measurement. See Details.
#' @param data optional data frame containing the data used to fit
#'   \code{object}. Used only when \code{object} is a fitted model object and
#'   the data cannot be recovered from it. Defaults to \code{NULL}.
#' @param maxiter positive integer giving the maximum number of iterations of
#'   the EM algorithm used to estimate the multivariate normal mean and
#'   covariance matrix. Defaults to \code{1000}.
#' @param tolerance a single positive finite number giving the relative
#'   convergence tolerance of the EM algorithm. Defaults to \code{1e-08}.
#' @param reference character string specifying the reference distribution:
#'   Little's exact F reference when the data are bivariate monotone and the
#'   chi-squared reference otherwise (\code{"auto"}), or always the chi-squared
#'   reference (\code{"asymptotic"}). Defaults to \code{"auto"}.
#'
#' @details
#' Let the rows be grouped into missing-data patterns. For pattern \eqn{r},
#' let \eqn{n_r} be its number of rows, \eqn{p_r} the number of observed
#' variables, and \eqn{\bar y_r}{ybar_r} their mean vector. With \eqn{\hat\mu}{mu} and
#' \eqn{\hat\Sigma}{Sigma} the maximum likelihood estimates of the common mean vector
#' and covariance matrix under a multivariate normal model, obtained by the EM
#' algorithm, and \eqn{\widetilde\Sigma = n\hat\Sigma/(n - 1)}{Sigma~ = n Sigma/(n - 1)}, Little's
#' statistic is
#' \deqn{d^2 = \sum_r n_r (\bar y_r - \hat\mu_r)^\top
#' \widetilde\Sigma_r^{-1}(\bar y_r - \hat\mu_r),}{d^2 = sum_r n_r (ybar_r - mu_r)' Sigma~_r^(-1) (ybar_r - mu_r),}
#' which under MCAR is asymptotically chi-squared with \eqn{\sum_r p_r - p}{sum_r p_r - p}
#' degrees of freedom, where \eqn{p} is the number of variables.
#'
#' When there are two variables, one fully observed, Little shows that
#' \eqn{d^2 = (n - 1)F/(n - 2 + F)} with \eqn{F \sim F_{1, n-2}}{F ~ F(1, n - 2)}. With
#' \code{reference = "auto"} this exact reference is used in that case,
#' provided the identity is reproduced numerically from the EM estimates;
#' otherwise a warning is issued and the chi-squared reference is used, which
#' usually indicates that \code{tolerance} is too loose.
#'
#' Matrix or data frame input must have at least two rows and two variables,
#' at least two observed values per variable, at least one jointly observed
#' value for every pair of variables, and no infinite values. Rows with no
#' observed values belong to no pattern and are removed with a warning. For a
#' fitted model object, the response is arranged with one row per cluster and
#' one column per \code{repeated} index; if \code{repeated} was not supplied
#' when fitting, the row order within each cluster defines the occasions.
#'
#' The test assumes a multivariate normal model and is most appropriate for
#' quantitative variables (Little, 1988), so a warning is issued for binary
#' variables. It is a diagnostic for the missingness mechanism rather than for
#' the fitted mean model (Hardin and Hilbe, 2013): a large p-value does not
#' establish MCAR. Without missing values the statistic is 0 and the p-value 1.
#'
#' @return
#' An object of class \code{c("geer_htest", "htest")}. The standard
#' components hold Little's statistic \eqn{d^2} (\code{statistic}) and the
#' degrees of freedom of its reference distribution (\code{parameter}):
#' \code{df} for the chi-squared reference, or \code{df1} and \code{df2} for
#' the exact F reference. In addition, the object contains:
#' \item{missing.patterns}{the number of distinct missing-data patterns.}
#' \item{n}{the number of independent rows or clusters tested.}
#' \item{variables}{the number of variables or repeated measurements tested.}
#' \item{iterations}{the number of EM iterations used.}
#' \item{converged}{whether the EM algorithm converged.}
#' \item{mean}{the maximum-likelihood estimate of the common mean vector.}
#' \item{covariance}{the degrees-of-freedom-corrected covariance matrix used
#' in Little's statistic.}
#' \item{ml.covariance}{the uncorrected multivariate-normal maximum-likelihood
#' covariance estimate.}
#' \item{reference}{the reference distribution used for the reported p-value.}
#' \item{asymptotic.df}{the degrees of freedom for the large-sample
#' chi-squared reference.}
#' \item{asymptotic.p.value}{the large-sample chi-squared p-value.}
#' \item{exact.p.value}{the exact bivariate monotone p-value when available,
#' otherwise \code{NULL}.}
#' \item{exact.f.statistic}{the corresponding bivariate ANOVA F statistic
#' when available, otherwise \code{NULL}.}
#'
#' @references
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' Little, R.J.A. (1988) A test of missing completely at random for
#' multivariate data with missing values. \emph{Journal of the American
#' Statistical Association}, \bold{83}, 1198--1202.
#'
#' @seealso \code{\link{mcar_logistic_test}},
#'   \code{\link{mcar_homoscedasticity_test}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' x <- data.frame(
#'   visit1 = c(2.1, 4.3, 3.2, 5.8, 7.1, 6.4, 8.2, 9.5, 10.1, 11.3),
#'   visit2 = c(1.5, 3.8, 2.9, 6.2, NA, 5.1, 8.9, NA, 9.4, 12.0),
#'   visit3 = c(3.0, 4.9, NA, 7.1, 6.3, 7.8, NA, 10.2, 11.5, 13.1)
#' )
#' mcar_little_test(x)
#'
#' @export
mcar_little_test <- function(object, data = NULL, maxiter = 1000L,
                             tolerance = 1e-8,
                             reference = c("auto", "asymptotic")) {
  data_name <- deparse1(substitute(object))
  maxiter <- check_integer_at_least(maxiter, "maxiter")
  if (!is_positive_scalar(tolerance)) {
    stop("'tolerance' must be a positive finite number", call. = FALSE)
  }
  reference <- match_mcar_reference(reference[1L])

  if (inherits(object, "geer")) {
    object <- check_geer_object(object)
    x <- extract_mcar_matrix(object, data = data)
    data_name <- "repeated response from fitted geer object"
  } else {
    if (!is.null(data)) {
      stop(
        "'data' is only used when 'object' is a fitted 'geer' object",
        call. = FALSE
      )
    }
    x <- object
  }

  x <- check_mcar_little_matrix(x)
  if (anyNA(x) && is_mcar_little_binary(x)) {
    warning(
      paste0(
        "Little (1988) notes that this MCAR test is most appropriate for ",
        "quantitative variables; consider contingency-table methods for ",
        "categorical variables"
      ),
      call. = FALSE
    )
  }

  result <- compute_mcar_little_statistic(
    x,
    maxiter = maxiter,
    tolerance = tolerance
  )
  exact <- result$bivariate_exact
  use_exact <- identical(reference, "auto") && !is.null(exact)

  if (use_exact) {
    statistic <- c("d^2" = result$statistic)
    parameter <- c(df1 = exact$df1, df2 = exact$df2)
    p_value <- exact$p_value
    method <- "Little's MCAR test (exact bivariate monotone reference)"
    reference_used <- "exact bivariate F"
  } else {
    statistic <- c("X-squared" = result$statistic)
    parameter <- c(df = result$df)
    p_value <- result$p_value
    method <- "Little's MCAR test"
    reference_used <- "asymptotic chi-squared"
  }

  structure(
    list(
      statistic = statistic,
      parameter = parameter,
      p.value = p_value,
      method = method,
      data.name = data_name,
      alternative = "data are not missing completely at random",
      missing.patterns = result$missing_patterns,
      n = nrow(x),
      variables = ncol(x),
      iterations = result$iterations,
      converged = result$converged,
      mean = result$mu,
      covariance = result$sigma,
      ml.covariance = result$sigma_ml,
      reference = reference_used,
      asymptotic.df = result$df,
      asymptotic.p.value = result$p_value,
      exact.p.value = if (is.null(exact)) NULL else exact$p_value,
      exact.f.statistic = if (is.null(exact)) NULL else exact$statistic
    ),
    class = c("geer_htest", "htest")
  )
}
