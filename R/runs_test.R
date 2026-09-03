#' @title
#' Wald-Wolfowitz Runs Test for geer Residuals
#'
#' @description
#' Performs the nonparametric Wald-Wolfowitz runs test proposed by Chang
#' (2000) for detecting non-random patterns in residuals from a fitted
#' generalized estimating equation model.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param type character string specifying the residual type. Options are
#'   \code{"pearson"}, \code{"deviance"}, and \code{"working"}. Defaults to
#'   \code{"pearson"}, matching the illustration in Hardin and Hilbe (2013).
#' @param order_by specifies the ordering of the residual sequence. The default,
#'   \code{"natural"}, orders observations by cluster identifier and then by
#'   the within-cluster \code{repeated} index. \code{"fitted"} orders by fitted
#'   values. A single model-matrix column name can be supplied to order by that
#'   covariate, or a numeric vector with one value per fitted observation can
#'   be supplied directly. Numeric vectors must already be aligned with the
#'   observation order stored in \code{object}. Ordering values must be finite
#'   and non-missing. Ties are resolved using the natural cluster/repeated
#'   order.
#'
#' @details
#' Let \eqn{n_p} and \eqn{n_n} denote the numbers of positive and negative
#' residuals, respectively, and let \eqn{T} be the observed number of runs in
#' their sign sequence. Under the null hypothesis that residual signs occur in
#' random order,
#' \deqn{E(T) = \frac{2n_p n_n}{n_p+n_n} + 1}
#' and
#' \deqn{V(T) =
#' \frac{2n_p n_n(2n_p n_n-n_p-n_n)}
#' {(n_p+n_n)^2(n_p+n_n-1)}.}
#' The standardized statistic
#' \deqn{Z = \frac{T-E(T)}{\sqrt{V(T)}}}
#' is compared with the standard normal distribution using a two-sided
#' p-value. No continuity correction is applied, following Chang (2000) and
#' Hardin and Hilbe (2013).
#'
#' Residuals equal to zero are omitted from the sign sequence because the test
#' is defined in terms of positive and negative residuals. Their number is
#' returned in the \code{zero} component. A warning is issued when either
#' \eqn{n_p <= 15} or \eqn{n_n <= 15}, since Chang (2000) recommends the
#' normal approximation for larger samples.
#'
#' The ordering is part of the hypothesis being assessed. Natural ordering can
#' diagnose remaining longitudinal or within-cluster structure. Ordering by a
#' continuous covariate can help assess its functional form, while ordering by
#' fitted values provides a broader model-adequacy diagnostic.
#'
#' @return
#' An object of class \code{"htest"}. In addition to the standard
#' \code{statistic}, \code{p.value}, \code{method}, \code{data.name}, and
#' \code{alternative} components, the object contains:
#' \item{runs}{the observed number of runs, \eqn{T}.}
#' \item{expected_runs}{the null expectation \eqn{E(T)}.}
#' \item{variance_runs}{the null variance \eqn{V(T)}.}
#' \item{positive}{the number of positive residuals, \eqn{n_p}.}
#' \item{negative}{the number of negative residuals, \eqn{n_n}.}
#' \item{zero}{the number of zero residuals omitted from the test.}
#' \item{nonzero}{the number of residuals used in the test.}
#' \item{residual_type}{the residual type used.}
#' \item{order_by}{a description of the ordering used.}
#'
#' @references
#' Chang, Y.-C. (2000) Residuals analysis of the generalized linear models for
#' longitudinal data. \emph{Statistics in Medicine}, \bold{19}, 1277--1293.
#' \doi{10.1002/(SICI)1097-0258(20000530)19:10<1277::AID-SIM494>3.0.CO;2-S}
#'
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' @seealso \code{\link{residuals.geer}}, \code{\link{fitted.geer}},
#'   \code{\link{predict.geer}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   seizures ~ treatment + lnbaseline + lnage,
#'   data = epilepsy,
#'   id = id,
#'   repeated = visit,
#'   family = poisson(link = "log"),
#'   corstr = "exchangeable"
#' )
#'
#' runs_test(fit)
#' runs_test(fit, type = "deviance")
#' runs_test(fit, order_by = "fitted")
#' runs_test(fit, order_by = "lnage")
#'
#' @export
runs_test <- function(object,
                      type = c("pearson", "deviance", "working"),
                      order_by = "natural") {
  object <- check_geer_object(object)
  type <- match.arg(type)
  ordering <- resolve_runs_order(object, order_by)
  residual_values <- stats::residuals(object, type = type)[ordering$index]
  result <- compute_runs_statistics(residual_values)

  if (result$positive <= 15L || result$negative <= 15L) {
    warning(
      "the normal approximation may be unreliable because there are 15 or fewer residuals of at least one sign",
      call. = FALSE
    )
  }

  structure(
    list(
      statistic = c(Z = result$statistic),
      p.value = result$p_value,
      method = "Wald-Wolfowitz runs test for geer residuals",
      data.name = sprintf("%s residuals ordered by %s", type, ordering$label),
      alternative = "two.sided",
      runs = result$runs,
      expected_runs = result$expected_runs,
      variance_runs = result$variance_runs,
      positive = result$positive,
      negative = result$negative,
      zero = result$zero,
      nonzero = result$nonzero,
      residual_type = type,
      order_by = ordering$label
    ),
    class = "htest"
  )
}
