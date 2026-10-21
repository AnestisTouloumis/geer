#' @title
#' Residual Goodness-of-Fit Tests for Binary geer Fits
#'
#' @description
#' Tests the fit of the marginal mean model of a binary response with the
#' Pearson chi-squared statistic or the unweighted sum of squared residuals,
#' standardized by their approximate means and variances (Pan, 2002).
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param statistic character string specifying the statistic: the Pearson
#'   chi-squared statistic (\code{"pearson"}) or the unweighted sum of squared
#'   residuals (\code{"uss"}). Defaults to \code{"pearson"}.
#' @param covariance character string specifying the estimator of the
#'   covariance matrix of the responses used in the variance of the
#'   statistic: the pooled unstructured estimator (\code{"unstructured"}) or
#'   the empirical estimator (\code{"empirical"}). See Details. Defaults to
#'   \code{"unstructured"}.
#'
#' @details
#' With fitted probabilities \eqn{\hat\mu_{ij}}{mu_ij}, the statistics are
#' \deqn{G = \sum_{i,j} \frac{(y_{ij} - \hat\mu_{ij})^2}{\hat\mu_{ij}(1 -
#' \hat\mu_{ij})}, \qquad U = \sum_{i,j} (y_{ij} - \hat\mu_{ij})^2,}{G =
#' sum_ij (y_ij - mu_ij)^2/(mu_ij (1 - mu_ij)),   U = sum_ij (y_ij - mu_ij)^2,}
#' with approximate means \eqn{N}, the number of observations, and
#' \eqn{\sum_{i,j} \hat\mu_{ij}(1 - \hat\mu_{ij})}{sum_ij mu_ij (1 - mu_ij)}.
#' Their approximate variances are \eqn{c^\top (I - H)\,\mathrm{cov}(Y)\,(I -
#' H)^\top c}{c' (I - H) cov(Y) (I - H)' c}, where \eqn{c} has elements
#' \eqn{(1 - 2\hat\mu_{ij})/(\hat\mu_{ij}(1 - \hat\mu_{ij}))}{(1 - 2 mu_ij)/
#' (mu_ij (1 - mu_ij))} for \eqn{G} and \eqn{1 - 2\hat\mu_{ij}}{1 - 2 mu_ij}
#' for \eqn{U}, and \eqn{H = D (\sum_i D_i^\top V_i^{-1} D_i)^{-1} D^\top
#' V^{-1}}{H = D (sum_i D_i' V_i^(-1) D_i)^(-1) D' V^(-1)} is the leverage
#' matrix of \code{\link{hatvalues.geer}}. The standardized statistic is
#' compared with the standard normal distribution, and the p-value is
#' two-sided, which reproduces the p-values reported by Pan (2002) and
#' Hosmer et al. (1997).
#'
#' With \code{covariance = "unstructured"}, \eqn{\mathrm{cov}(Y_i)}{cov(Y_i)}
#' is estimated by \eqn{A_i^{1/2} \hat R A_i^{1/2}}{A_i^(1/2) R A_i^(1/2)},
#' where \eqn{A_i} is the diagonal matrix of the
#' \eqn{\hat\mu_{ij}(1 - \hat\mu_{ij})}{mu_ij (1 - mu_ij)} and \eqn{\hat R}{R}
#' the pooled correlation matrix of the standardized residuals. Each element
#' of \eqn{\hat R}{R} is the average of the products of the standardized
#' residuals at the two \code{repeated} positions over the clusters that
#' observe both. Unlike the pooled estimators of \code{\link{vcov.geer}},
#' such as \code{cov_type = "pan"}, which require a balanced design,
#' \eqn{\hat R}{R} is therefore available for unbalanced designs. Pan (2002)
#' studied balanced designs, and with unbalanced designs the positions should
#' be given by \code{repeated}. With \code{"empirical"}, the covariance is
#' estimated by \eqn{(y_i - \hat\mu_i)(y_i - \hat\mu_i)^\top}{(y_i - mu_i)
#' (y_i - mu_i)'}. Pan (2002) recommends the unstructured estimator. Pan
#' (2002) and Evans and Li (2005) label the two versions of each statistic in
#' opposite ways.
#'
#' Pan (2002) derives the variances for the logit link and recommends a
#' working independence fit, under which the approximation of the residuals
#' is most accurate. Here the working covariance matrices \eqn{V_i} of
#' \code{object} are used, and \eqn{D_i = \partial\mu_i/\partial\beta^\top}{D_i
#' = d mu_i/d beta'} for any link, so \code{object} can be a
#' \code{\link{geewa}} fit of a binomial family or a
#' \code{\link{geewa_binary}} fit. The response must be binary with unit
#' prior weights. The tests are intended for data with many covariate
#' patterns, such as when the model contains a continuous covariate. With few
#' covariate patterns, the variance of the statistic is close to zero, the
#' Pearson and unweighted sum of squares statistics can give the same
#' p-value, and the tests are not informative.
#'
#' The tests were developed for \code{method = "gee"}. For the other
#' estimation methods, the statistic is computed from the fitted values of
#' \code{object} and a warning is issued, since its approximate mean and
#' variance have not been studied for these estimators.
#'
#' @return
#' An object of class \code{c("geer_htest", "htest")}. The standard
#' components hold the standardized statistic (\code{statistic}), the
#' statistic (\code{estimate}), and its approximate mean
#' (\code{null.value}). In addition, the object contains:
#' \item{mean}{the approximate mean of the statistic.}
#' \item{variance}{the approximate variance of the statistic.}
#' \item{covariance}{the estimator of the covariance matrix of the responses
#'   used.}
#'
#' @references
#' Evans, S. and Li, L. (2005) A comparison of goodness of fit tests for the
#' logistic GEE model. \emph{Statistics in Medicine}, \bold{24}, 1245--1261.
#'
#' Hosmer, D.W., Hosmer, T., le Cessie, S. and Lemeshow, S. (1997) A
#' comparison of goodness-of-fit tests for the logistic regression model.
#' \emph{Statistics in Medicine}, \bold{16}, 965--980.
#'
#' Pan, W. (2002) Goodness-of-fit tests for GEE with correlated binary data.
#' \emph{Scandinavian Journal of Statistics}, \bold{29}, 101--110.
#'
#' @seealso \code{\link{gof_score_test}}, \code{\link{gof_tests}},
#'   \code{\link{hatvalues.geer}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", ]
#' fit <- geewa_binary(
#'   formula = status ~ baseline + treatment + gender + visit + age,
#'   link = "logit",
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   orstr = "independence"
#' )
#' gof_residual_test(fit)
#' gof_residual_test(fit, statistic = "uss", covariance = "empirical")
#'
#' @export
gof_residual_test <- function(object,
                              statistic = c("pearson", "uss"),
                              covariance = c("unstructured", "empirical")) {
  context <- "gof_residual_test"
  data_name <- deparse1(substitute(object))
  object <- check_geer_object(object)
  statistic <- match_gof_statistic(statistic[1L])
  covariance <- match_gof_covariance(covariance[1L])
  check_gof_fit(object, context)
  result <- compute_gof_residual_test(
    compute_gof_residual_parts(object, covariance),
    statistic,
    context
  )
  label <- if (identical(statistic, "pearson")) {
    "Pearson statistic"
  } else {
    "unweighted sum of squares"
  }
  structure(
    list(
      statistic = c(Z = result$z),
      estimate = stats::setNames(result$value, label),
      null.value = stats::setNames(result$mean, label),
      p.value = result$p_value,
      alternative = "two.sided",
      method = sprintf(
        "Pan goodness-of-fit test (%s) with the %s covariance matrix",
        label,
        covariance
      ),
      data.name = data_name,
      mean = result$mean,
      variance = result$variance,
      covariance = covariance
    ),
    class = c("geer_htest", "htest")
  )
}
