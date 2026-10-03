#' @title
#' Little's Test for Missing Completely at Random Data
#'
#' @description
#' Performs Little's (1988) test of the null hypothesis that
#' missing values are missing completely at random (MCAR). For a fitted
#' \code{geer} object, the repeated response is reconstructed in wide form
#' from the original data before rows with missing responses were omitted by
#' model fitting.
#'
#' @param object a fitted object of class \code{"geer"}, a numeric matrix, or
#'   a numeric data frame. For matrix or data-frame input, rows are independent
#'   units and columns are variables or repeated measurements. The supplied
#'   data must contain at least two rows and two variables, with at least two
#'   observed values per variable and at least one jointly observed value for
#'   every pair of variables. Infinite values are not allowed.
#' @param data optional original data used to fit \code{object}. This is only
#'   used when \code{object} is a \code{geer} fit and is useful when the
#'   original data cannot be recovered from the fitted object. Defaults to
#'   \code{NULL}.
#' @param maxit positive integer giving the maximum number of EM iterations
#'   used to obtain the multivariate-normal maximum-likelihood estimates.
#'   Defaults to 1000.
#' @param tol a single positive finite number specifying the convergence
#'   tolerance for the EM algorithm. Defaults to \code{1e-08}.
#' @param reference reference distribution used for the p-value. \code{"auto"}
#'   uses Little's exact bivariate monotone result when it applies and
#'   otherwise falls back, without warning, to the asymptotic chi-squared
#'   distribution. \code{"asymptotic"} always uses the chi-squared reference.
#'   Defaults to \code{"auto"}.
#'
#' @details
#' Let the rows be divided into missing-data patterns. For pattern \eqn{r},
#' let \eqn{n_r} be its number of rows and \eqn{\bar y_r} the mean vector
#' of the variables observed in that pattern. Let \eqn{\hat\mu} and
#' \eqn{\hat\Sigma} be the common mean vector and covariance matrix estimated
#' by maximum likelihood under an ignorable multivariate-normal missing-data
#' model, and define Little's degrees-of-freedom correction
#' \eqn{\widetilde\Sigma = n\hat\Sigma/(n-1)}. Little's statistic is
#' \deqn{d^2 = \sum_r n_r (\bar y_r-\hat\mu_r)^T
#' \widetilde\Sigma_r^{-1}(\bar y_r-\hat\mu_r).}
#' Under MCAR, the statistic is asymptotically chi-squared with
#' \deqn{\mathrm{df} = \sum_r p_r - p,}
#' where \eqn{p_r} is the number of variables observed in pattern \eqn{r}
#' and \eqn{p} is the total number of variables.
#'
#' The common mean and covariance are estimated internally with an EM
#' algorithm, so no additional package is required. Following Little (1988),
#' the multivariate-normal maximum-likelihood covariance estimate
#' \eqn{\hat\Sigma} is multiplied by \eqn{n/(n-1)} before its pattern-specific
#' submatrices are used in the test statistic.
#'
#' For two variables with one variable observed for every row and missingness
#' confined to the other variable, Little shows that the small-sample null
#' distribution can be obtained from an ordinary two-group ANOVA (equivalently,
#' a pooled two-sample t test):
#' \deqn{d^2 = \frac{(n-1)F}{n-2+F}, \qquad F \sim F_{1,n-2}.}
#' With \code{reference = "auto"}, this exact normal-theory reference is used
#' automatically when the data fall in that special bivariate monotone case
#' and the corresponding Little-statistic identity is reproduced to numerical
#' tolerance. If the case is detected but the identity is not reproduced, a
#' warning is issued and the asymptotic chi-squared reference is used instead.
#' The identity is checked against the statistic obtained from the EM
#' estimates, so it is sensitive to \code{tol}: at the default the two agree
#' to well within the comparison tolerance, whereas a substantially looser
#' \code{tol} can prevent the identity from being reproduced and therefore
#' cost the exact reference. A warning here indicates insufficient convergence
#' rather than an inapplicable reference.
#' For all other missing-data patterns the large-sample chi-squared reference
#' is used silently, with no warning. Little also derives a more general
#' small-sample distribution for monotone patterns as a sum of transformed
#' independent F variables; that nonstandard reference distribution is not
#' evaluated here.
#'
#' Rows with no observed values belong to no missing-data pattern in Little's
#' construction. They are removed, with a warning, before the statistic and
#' the sample size \eqn{n} are computed.
#'
#' For a \code{geer} fit, only the repeated response is tested. The function
#' reconstructs a cluster-by-repeated-measure matrix from \code{id} and
#' \code{repeated}. If \code{repeated} was omitted during fitting, the
#' within-cluster row order in the original data defines the repeated-measure
#' positions. Consequently, an entirely absent row cannot be distinguished
#' from a measurement occasion that never existed unless \code{repeated} was
#' supplied explicitly.
#'
#' Little's test is based on a multivariate-normal working model for the
#' variables being tested. Little (1988) notes that the test is most appropriate
#' for quantitative variables and recommends contingency-table methods for
#' categorical variables. When missing values are present, the function warns
#' if a binary variable is detected. The paper also reports that the asymptotic
#' chi-squared test can
#' be conservative in small samples.
#'
#' As emphasized by Hardin and Hilbe (2013), this is a diagnostic for the
#' missingness mechanism rather than a test of the fitted GEE mean model itself.
#' A small p-value provides evidence against MCAR; a large p-value does not
#' establish that MCAR holds. If the supplied data contain no missing values,
#' the function returns a test statistic of 0 with p-value 1.
#'
#' @return
#' An object of class \code{"htest"}. In addition to the standard components,
#' the object contains:
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
#' @seealso \code{\link{runs_test}}, \code{\link{residuals.geer}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
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
mcar_little_test <- function(object, data = NULL, maxit = 1000L, tol = 1e-8,
                             reference = c("auto", "asymptotic")) {
  data_name <- deparse1(substitute(object))
  maxit <- check_integer_at_least(maxit, "maxit")
  if (!is_positive_scalar(tol)) {
    stop("'tol' must be a single positive finite number", call. = FALSE)
  }
  reference <- match.arg(reference)

  if (inherits(object, "geer")) {
    object <- check_geer_object(object)
    x <- extract_geer_mcar_matrix(object, data = data)
    data_name <- "repeated response from fitted geer object"
  } else {
    if (!is.null(data)) {
      stop("'data' is only used when 'object' is a fitted 'geer' object", call. = FALSE)
    }
    x <- object
  }

  x <- check_mcar_little_matrix(x)
  if (anyNA(x) && has_mcar_little_binary_variable(x)) {
    warning(
      paste0(
        "Little (1988) notes that this MCAR test is most appropriate for ",
        "quantitative variables; consider contingency-table methods for ",
        "categorical variables"
      ),
      call. = FALSE
    )
  }

  result <- compute_mcar_little_statistic(x, maxit = maxit, tol = tol)
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
    class = "htest"
  )
}


