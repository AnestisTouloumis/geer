#' @title
#' Wald-Wolfowitz Runs Test for geer Residuals
#'
#' @description
#' Performs the nonparametric Wald-Wolfowitz runs test proposed by Chang
#' (2000) for detecting nonrandom patterns in residuals from a fitted
#' generalized estimating equation model.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param order_by character string or numeric vector specifying the order of
#'   the residuals: by cluster and \code{repeated} index (\code{"natural"}), by
#'   fitted value (\code{"fitted"}), by a covariate named in the model matrix or
#'   in the data used to fit the model, or by a numeric vector with one finite
#'   value per observation, in the order stored in \code{object}. Ties are
#'   broken by the natural order. See Details. Defaults to \code{"natural"}.
#' @param alternative character string specifying the alternative hypothesis:
#'   \code{"two.sided"}, \code{"less"} (fewer runs than expected), or
#'   \code{"greater"} (more runs than expected). Defaults to \code{"two.sided"}.
#'
#' @details
#' Let \eqn{n_p} and \eqn{n_n} be the numbers of positive and negative
#' residuals and \eqn{T} the number of runs in their sign sequence. Under the
#' null hypothesis that the signs occur in random order,
#' \deqn{E(T) = \frac{2n_p n_n}{n_p + n_n} + 1, \qquad
#' V(T) = \frac{2n_p n_n(2n_p n_n - n_p - n_n)}{(n_p + n_n)^2(n_p + n_n - 1)},}{E(T) = 2 n_p n_n/(n_p + n_n) + 1,   V(T) = 2 n_p n_n (2 n_p n_n - n_p - n_n)/((n_p + n_n)^2 (n_p + n_n - 1)),}
#' and \eqn{Z = (T - E(T))/\sqrt{V(T)}}{Z = (T - E(T))/sqrt(V(T))} is compared with the standard normal
#' distribution without continuity correction, as in Chang (2000) and Hardin
#' and Hilbe (2013). Chang (2000) recommends the approximation only when
#' \eqn{n_p > 15} and \eqn{n_n > 15}; otherwise a warning is issued.
#'
#' \code{alternative = "less"} tests against too few runs, as produced by
#' clustered signs, an unmodeled trend, or omitted structure, and
#' \code{"greater"} against too many runs, that is, alternating signs. Chang
#' (2000) rejects for extreme \eqn{Z} in either direction, which corresponds
#' to the two-sided default. Hardin and Hilbe (2013) report one-sided p-values
#' in the direction of the observed \eqn{Z}.
#'
#' Only the signs of the residuals are used. The working, Pearson, and
#' deviance residuals have the same signs, so the test does not depend on the
#' residual type. Zero residuals are omitted and counted in the \code{zero}
#' component. Prior weights are ignored: each observation contributes one sign.
#'
#' A character \code{order_by} that is not \code{"natural"}, \code{"fitted"},
#' or a model matrix column is looked up in the data used to fit the model,
#' so a covariate omitted from the model can be used. The model frame is
#' rebuilt from the call and checked against the stored clusters and
#' occasions; an error is signaled if the variable adds missing values or the
#' check fails, in which case an aligned numeric vector can be supplied
#' instead. Factor, character, and logical variables are ordered by their
#' codes, and date-time variables by their numeric values.
#'
#' @section Choice of ordering:
#' The null distribution does not depend on the ordering; the ordering
#' determines the alternative. Ordering by a variable gives power against
#' signs that depend on that variable:
#' \describe{
#'   \item{\code{"natural"}}{signs that are grouped within clusters or
#'     drift over time;}
#'   \item{\code{"fitted"}}{signs that depend on the fitted mean, suggesting
#'     an inadequate link or variance function;}
#'   \item{a model covariate}{a misspecified functional form of that
#'     covariate;}
#'   \item{another covariate}{a wrongly omitted covariate;}
#'   \item{a numeric vector}{whatever the vector indexes, such as batch or
#'     center.}
#' }
#' Chang (2000) and Hardin and Hilbe (2013, Section 4.2.1) use the natural
#' ordering, and Hardin and Hilbe also order by a covariate; ordering by
#' fitted values or by a supplied vector is an extension.
#'
#' Under the natural ordering, too few runs can arise from positive
#' within-cluster association as well as from a misspecified mean, so a
#' rejection does not identify the cause. Orderings that vary within clusters
#' interleave the clusters and make a rejection more readily attributable to
#' the mean structure. A cluster-level covariate, such as treatment arm, keeps
#' the clusters together and behaves like the natural ordering.
#'
#' Ties in the ordering are broken by the natural order, so a heavily tied
#' ordering, such as fitted values from few covariate patterns or a factor
#' with few levels, approaches the natural ordering.
#'
#' Hardin and Hilbe (2013, Section 4.2.1) print \eqn{E(T) = 36.32} for their
#' second example under the natural ordering; the correct value is
#' \eqn{2(22)(78)/100 + 1 = 35.32}, which gives \eqn{Z = 0.494} instead of
#' \eqn{Z = 0.200}.
#'
#' @return
#' An object of class \code{c("geer_runs_test", "geer_htest", "htest")}, which
#' can be plotted with \code{\link{plot.geer_runs_test}}. The standard
#' components hold \eqn{Z} (\code{statistic}), \eqn{n_p} and \eqn{n_n}
#' (\code{parameter}), \eqn{T} (\code{estimate}), and \eqn{E(T)}
#' (\code{null.value}). In addition, the object contains:
#' \item{runs}{the observed number of runs, \eqn{T}.}
#' \item{expected_runs}{the null expectation \eqn{E(T)}.}
#' \item{variance_runs}{the null variance \eqn{V(T)}.}
#' \item{positive}{the number of positive residuals, \eqn{n_p}.}
#' \item{negative}{the number of negative residuals, \eqn{n_n}.}
#' \item{zero}{the number of zero residuals omitted from the test.}
#' \item{nonzero}{the number of residuals used in the test.}
#' \item{order_by}{a description of the ordering used.}
#' \item{natural_order}{a logical value indicating whether the natural
#'   cluster/repeated ordering was used.}
#' \item{signs}{the sign sequence whose runs were counted, in the order
#'   tested and with zero residuals removed.}
#' \item{cluster}{the cluster identifier of each element of \code{signs}.}
#'
#' @references
#' Chang, Y.-C. (2000) Residuals analysis of the generalized linear models for
#' longitudinal data. \emph{Statistics in Medicine}, \bold{19}, 1277--1293.
#'
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' @seealso \code{\link{plot.geer_runs_test}}, \code{\link{residuals.geer}},
#'   \code{\link{fitted.geer}}, \code{\link{predict.geer}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   repeated = visit,
#'   corstr = "exchangeable"
#' )
#' runs_test(fit)
#' runs_test(fit, order_by = "fitted")
#' runs_test(fit, order_by = "lnage")
#' runs_test(fit, alternative = "less")
#'
#' @export
runs_test <- function(object,
                      order_by = "natural",
                      alternative = c("two.sided", "less", "greater")) {
  caller_env <- parent.frame()
  object <- check_geer_object(object)
  alternative <- match_alternative(alternative[1L])
  ordering <- get_runs_order(object, order_by, env = caller_env)
  residual_values <-
    stats::residuals(object, type = "working")[ordering$index]
  result <- compute_runs_statistics(residual_values, alternative)
  cluster <- object$id[ordering$index][result$retained]
  if (result$positive <= 15L || result$negative <= 15L) {
    warning(
      paste0(
        "the normal approximation may be unreliable because there are 15 or ",
        "fewer residuals of at least one sign"
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      statistic = c(Z = result$statistic),
      parameter = c(n_p = result$positive, n_n = result$negative),
      p.value = result$p_value,
      estimate = c("number of runs" = as.numeric(result$runs)),
      null.value = c("number of runs" = result$expected_runs),
      method = "Wald-Wolfowitz runs test for geer residuals",
      data.name = sprintf("residuals ordered by %s", ordering$label),
      alternative = alternative,
      runs = result$runs,
      expected_runs = result$expected_runs,
      variance_runs = result$variance_runs,
      positive = result$positive,
      negative = result$negative,
      zero = result$zero,
      nonzero = result$nonzero,
      order_by = ordering$label,
      natural_order = isTRUE(ordering$natural),
      signs = result$signs,
      cluster = cluster
    ),
    class = c("geer_runs_test", "geer_htest", "htest")
  )
}
