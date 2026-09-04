#' @title
#' Wald-Wolfowitz Runs Test for geer Residuals
#'
#' @description
#' Performs the nonparametric Wald-Wolfowitz runs test proposed by Chang
#' (2000) for detecting non-random patterns in residuals from a fitted
#' generalized estimating equation model.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param order_by specifies the ordering of the residual sequence. The default,
#'   \code{"natural"}, orders observations by cluster identifier and then by
#'   the within-cluster \code{repeated} index. \code{"fitted"} orders by fitted
#'   values. A single model-matrix column name can be supplied to order by that
#'   covariate, or a numeric vector with one value per fitted observation can
#'   be supplied directly. Numeric vectors must already be aligned with the
#'   observation order stored in \code{object}. A name that is not a
#'   model-matrix column is looked up among the variables of the data used to
#'   fit the model, so a covariate omitted from the model can be used; see
#'   Details. Ordering values must be finite
#'   and non-missing. Ties are resolved using the natural cluster/repeated
#'   order; see Details for the consequences when the ordering has few
#'   distinct values.
#' @param alternative character string specifying the alternative hypothesis,
#'   expressed in terms of the observed number of runs \eqn{T} relative to its
#'   null expectation. Options are \code{"two.sided"}, \code{"less"} (fewer
#'   runs than expected) and \code{"greater"} (more runs than expected).
#'   Defaults to \code{"two.sided"}.
#' @param nperm \code{NULL} or a single positive whole number. The default,
#'   \code{NULL}, refers the standardized statistic to the normal
#'   approximation of Chang (2000). A number instead requests a Monte Carlo
#'   p-value obtained from that many within-cluster permutations of the
#'   residual signs; see Details.
#' @param exact \code{NULL} or a single logical value controlling whether the
#'   exact null distribution of \eqn{T} is used instead of the normal
#'   approximation. The default, \code{NULL}, uses the exact distribution
#'   whenever the normal approximation is not recommended, that is when
#'   \eqn{n_p <= 15} or \eqn{n_n <= 15}. Cannot be combined with
#'   \code{nperm}. It is ignored, without an error, when \code{nperm} is
#'   supplied: a permutation reference distribution takes precedence over both
#'   the exact and the approximate one.
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
#' is compared with the standard normal distribution. No continuity correction
#' is applied, following Chang (2000) and Hardin and Hilbe (2013).
#'
#' The \code{alternative} argument selects the tail of that distribution.
#' \code{"less"} tests against too few runs, the alternative implied by sign
#' clustering, an unmodeled trend, or omitted structure, and is the departure
#' emphasized by Chang (2000). \code{"greater"} tests against too many runs,
#' that is, systematic alternation of residual signs. Under the natural
#' ordering a small lower-tail p-value is ambiguous: too few runs arise both
#' from a misspecified mean structure and from genuine positive within-cluster
#' association, because the null distribution of \eqn{T} treats the signs as
#' exchangeable across observations. The default is therefore two-sided.
#'
#' Hardin and Hilbe (2013) report a one-sided p-value for their worked
#' example, so reproducing their Section 4.2.1 figure requires
#' \code{alternative = "greater"}; the default two-sided p-value is twice as
#' large.
#'
#' The test uses the residuals only through their signs, and therefore does
#' not depend on which residual type is extracted from the fitted model. The
#' working (raw), Pearson and deviance residuals returned by
#' \code{\link{residuals.geer}} are, respectively, \eqn{y_{ij}-\hat\mu_{ij}},
#' that difference divided by a strictly positive scale factor, and its signed
#' square root. All three therefore produce the same sign sequence, the same
#' \eqn{n_p}, \eqn{n_n} and \eqn{T}, and hence the same value of \eqn{Z} and
#' the same p-value, so no residual-type argument is offered. This invariance
#' is also why Chang (2000), who works with deviance residuals, and Hardin and
#' Hilbe (2013), who plot raw residuals, describe the same test.
#'
#' Residuals equal to zero are omitted from the sign sequence because the test
#' is defined in terms of positive and negative residuals. Their number is
#' returned in the \code{zero} component.
#'
#' Prior weights play no part in the test. Every retained observation
#' contributes exactly one sign, so a row of a two-column binomial response
#' representing many trials counts once, as does an observation with a large
#' analytic weight.
#'
#' Chang (2000) recommends the normal approximation only when \eqn{n_p > 15}
#' and \eqn{n_n > 15}. Outside that range \code{runs_test} evaluates the
#' exact null distribution of \eqn{T} instead, which for \eqn{T = 2k} and
#' \eqn{T = 2k+1} respectively is
#' \deqn{P(T = 2k) = \frac{2\binom{n_p-1}{k-1}\binom{n_n-1}{k-1}}
#' {\binom{n_p+n_n}{n_p}}, \qquad
#' P(T = 2k+1) = \frac{\binom{n_p-1}{k}\binom{n_n-1}{k-1} +
#' \binom{n_p-1}{k-1}\binom{n_n-1}{k}}{\binom{n_p+n_n}{n_p}}.}
#' Tail probabilities are summed over this distribution and, for a two-sided
#' test, the smaller tail is doubled and truncated at one. Setting
#' \code{exact = TRUE} or \code{exact = FALSE} forces or suppresses this
#' behavior; when the normal approximation is used with small sign counts a
#' warning is issued instead.
#'
#' Exact and approximate p-values do not agree closely even at large sign
#' counts, because no continuity correction is applied to the normal
#' approximation. The exact tail probabilities match the continuity-corrected
#' approximation instead, and the difference between the two is of that order:
#' with \eqn{n_p = n_n = 500} and \eqn{T = 520}, for instance, the exact
#' upper-tail probability is about \eqn{0.1209} against \eqn{0.1146} from the
#' uncorrected approximation. The exact value should be preferred when the two
#' are compared.
#'
#' The ordering is part of the hypothesis being assessed. Chang (2000)
#' considers only the natural ordering, in which each cluster's measurements
#' appear consecutively and in ascending time order, and uses it to diagnose
#' remaining longitudinal or within-cluster structure. Ordering by fitted
#' values or by a covariate is the amendment of Hardin and Hilbe (2013,
#' Section 4.2.1): ordering by a continuous covariate can help assess its
#' functional form, while ordering by fitted values provides a broader
#' model-adequacy diagnostic.
#'
#' Orderings other than the natural one are only well defined up to ties, and
#' ties are common rather than exceptional. A factor column takes as many
#' distinct values as it has levels, and fitted values take one distinct value
#' per covariate pattern, so a design with few covariate patterns yields few
#' distinct ordering values; the worked example of Hardin and Hilbe (2013) has
#' three distinct fitted values across eighty observations. Within a group of
#' tied ordering values the sequence, and therefore \eqn{T}, is arbitrary.
#' Ties are broken here by the natural cluster/repeated order, which makes the
#' result reproducible but means that a heavily tied ordering reduces to the
#' natural ordering within each tie group. Such an ordering should be
#' interpreted accordingly, and a statistic computed from an ordering with few
#' distinct values carries correspondingly little information about the
#' quantity being ordered on.
#'
#' A character \code{order_by} that names neither \code{"natural"},
#' \code{"fitted"}, nor a column of the model matrix is resolved against the
#' data used to fit the model, which allows the residuals to be ordered on a
#' covariate that was left out of the model. Such a variable cannot simply be
#' taken from the stored data frame, because rows are dropped by
#' \code{na.action} and by \code{subset} and the retained rows are then
#' sorted by cluster and time before being stored. The model frame is instead
#' rebuilt from the recorded call with the variable added to the formula, and
#' the reconstruction is checked against the stored cluster and time indices.
#' An error is signalled, rather than a possibly misaligned ordering returned,
#' if the variable introduces additional missing values or if the check fails;
#' an aligned numeric vector can always be supplied directly instead. Factor
#' and character variables are ordered by their factor codes.
#'
#' Chang (2000) further notes that for models with Poisson responses the
#' residual sequence should follow the order of measurement rather than the
#' follow-up time, which the default natural ordering provides.
#'
#' \strong{Permutation p-values.} The null distribution of \eqn{T} given
#' above treats the residual signs as exchangeable across all observations,
#' that is, as if the observations were independent. Under the natural
#' ordering this is violated by the within-cluster association that the model
#' is fitted to accommodate: positively associated residuals produce runs of
#' a common sign, \eqn{T} is too small, and the test can reject even when the
#' mean structure is correct. Chang (2000) treats such a rejection as grounds
#' for abandoning the marginal model in favor of a random-effects model, which
#' changes the estimand from marginal to conditional.
#'
#' Supplying \code{nperm} replaces the normal approximation by a Monte Carlo
#' test that does not require exchangeability across clusters. The residual
#' signs are permuted independently within each cluster, the sequence is
#' re-read in the requested ordering, and \eqn{T} is recomputed; the reported
#' p-value is \eqn{(1 + m) / (\code{nperm} + 1)}, where \eqn{m} counts the
#' permutations at least as extreme as the observed \eqn{T} in the direction
#' given by \code{alternative}. Two-sided p-values are twice the smaller tail
#' probability, truncated at one. This is a deliberate extension of Chang
#' (2000), who considers only the normal approximation, and the two p-values
#' will not agree.
#'
#' The permuted null is that residual signs are exchangeable \emph{within}
#' each cluster. A cluster whose residuals all share one sign contributes the
#' same number of runs under every permutation, so cluster-level lopsidedness
#' can no longer produce a rejection on its own, while an ordering-related
#' pattern within clusters still can. Two consequences follow. If most
#' clusters are internally constant in sign the permutation distribution is
#' nearly degenerate and the test has little power; a warning is issued if it
#' is exactly degenerate. And serial dependence is still confounded with mean
#' misspecification, since autoregressive residuals are not exchangeable
#' within a cluster either; the permutation test removes the between-cluster
#' component of the problem, not the within-cluster one.
#'
#' The reported statistic follows the null that was used. With \code{nperm}
#' the number of runs is standardized by the mean and standard deviation of
#' the permutation distribution rather than by \eqn{E(T)} and \eqn{V(T)},
#' and \code{null.value} reports the permutation mean, so that the printed
#' statistic and p-value refer to the same reference distribution. The
#' statistic is \code{NA} when the permutation distribution is degenerate.
#' The analytic moments remain available in the \code{expected_runs} and
#' \code{variance_runs} components whichever null is used.
#'
#' Permutation p-values are random. Use \code{\link{set.seed}} for
#' reproducible results, and note that the attainable resolution is
#' \eqn{1/(\code{nperm}+1)}. The small-sample warning about \eqn{n_p} and
#' \eqn{n_n} is not issued when \code{nperm} is supplied, since the normal
#' approximation is then not used.
#'
#' @return
#' An object of class \code{c("geer_runs_test", "htest")}, so that it prints
#' through \code{\link[stats]{print.htest}} and can be drawn with
#' \code{\link{plot.geer_runs_test}}. The standard components are used as
#' follows: \code{statistic} holds \eqn{Z}, \code{parameter} holds
#' \eqn{n_p} and \eqn{n_n}, \code{estimate} holds the observed number of
#' runs and \code{null.value} its null expectation, so that all of the
#' quantities entering the test are shown by the default print method. In
#' addition, the object contains:
#' \item{runs}{the observed number of runs, \eqn{T}.}
#' \item{expected_runs}{the exchangeable null expectation \eqn{E(T)}.}
#' \item{variance_runs}{the exchangeable null variance \eqn{V(T)}.}
#' \item{permuted_mean, permuted_sd}{the mean and standard deviation of the
#'   permutation distribution of \eqn{T}, or \code{NA} when no permutation
#'   test was performed.}
#' \item{positive}{the number of positive residuals, \eqn{n_p}.}
#' \item{negative}{the number of negative residuals, \eqn{n_n}.}
#' \item{zero}{the number of zero residuals omitted from the test.}
#' \item{nonzero}{the number of residuals used in the test.}
#' \item{nperm}{the number of within-cluster permutations used, or \code{NA}
#'   when no permutation test was performed.}
#' \item{exact}{a logical value indicating whether the exact null
#'   distribution of \eqn{T} was used.}
#' \item{natural_order}{a logical value indicating whether the natural
#'   cluster/repeated ordering was used.}
#' \item{signs}{the sign sequence whose runs were counted, in the order
#'   tested and with zero residuals removed.}
#' \item{cluster}{the cluster identifier of each element of \code{signs}.}
#' \item{order_by}{a description of the ordering used.}
#'
#' @references
#' Chang, Y.-C. (2000) Residuals analysis of the generalized linear models for
#' longitudinal data. \emph{Statistics in Medicine}, \bold{19}, 1277--1293.
#'
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' @seealso \code{\link{plot.geer_runs_test}},
#'   \code{\link{residuals.geer}}, \code{\link{fitted.geer}},
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
#' runs_test(fit, order_by = "fitted")
#' runs_test(fit, order_by = "lnage")
#' runs_test(fit, alternative = "less")
#'
#' set.seed(1)
#' runs_test(fit, nperm = 999)
#'
#' @export
runs_test <- function(object,
                      order_by = "natural",
                      alternative = c("two.sided", "less", "greater"),
                      nperm = NULL,
                      exact = NULL) {
  caller_env <- parent.frame()
  object <- check_geer_object(object)
  alternative <- match.arg(alternative)
  nperm <- normalize_runs_nperm(nperm)
  exact <- normalize_runs_exact(exact, nperm)
  ordering <- resolve_runs_order(object, order_by, env = caller_env)
  residual_values <-
    stats::residuals(object, type = "working")[ordering$index]
  result <- compute_runs_statistics(residual_values, alternative)
  retained <- sign(residual_values) != 0
  signs <- as.integer(sign(residual_values)[retained])
  cluster <- as.numeric(object$id[ordering$index])[retained]
  small_counts <- result$positive <= 15L || result$negative <= 15L
  use_exact <- if (is.null(exact)) is.null(nperm) && small_counts else exact
  statistic_value <- result$statistic
  null_runs <- result$expected_runs
  permuted_mean <- NA_real_
  permuted_sd <- NA_real_
  if (!is.null(nperm)) {
    permutation <- compute_runs_permutation(
      residual_values = residual_values,
      cluster = object$id[ordering$index],
      runs = result$runs,
      alternative = alternative,
      nperm = nperm
    )
    result$p_value <- permutation$p_value
    permuted_mean <- permutation$mean
    permuted_sd <- permutation$sd
    null_runs <- permuted_mean
    statistic_value <- if (is.finite(permuted_sd) && permuted_sd > 0) {
      (result$runs - permuted_mean) / permuted_sd
    } else {
      NA_real_
    }
    method <- sprintf(
      "Wald-Wolfowitz runs test for geer residuals (%d within-cluster permutations)",
      nperm
    )
  } else if (use_exact) {
    result$p_value <- compute_runs_exact_p(
      positive_no = result$positive,
      negative_no = result$negative,
      runs = result$runs,
      alternative = alternative
    )
    method <- "Wald-Wolfowitz runs test for geer residuals (exact null distribution)"
  } else {
    method <- "Wald-Wolfowitz runs test for geer residuals"
    if (small_counts) {
      warning(
        "the normal approximation may be unreliable because there are 15 or fewer residuals of at least one sign",
        call. = FALSE
      )
    }
  }

  structure(
    list(
      statistic = c(Z = statistic_value),
      parameter = c(n_p = result$positive, n_n = result$negative),
      p.value = result$p_value,
      estimate = c("number of runs" = as.numeric(result$runs)),
      null.value = c("number of runs" = null_runs),
      method = method,
      data.name = sprintf("residuals ordered by %s", ordering$label),
      alternative = alternative,
      runs = result$runs,
      expected_runs = result$expected_runs,
      variance_runs = result$variance_runs,
      permuted_mean = permuted_mean,
      permuted_sd = permuted_sd,
      positive = result$positive,
      negative = result$negative,
      zero = result$zero,
      nonzero = result$nonzero,
      nperm = if (is.null(nperm)) NA_integer_ else nperm,
      exact = isTRUE(use_exact),
      order_by = ordering$label,
      natural_order = isTRUE(ordering$natural),
      signs = signs,
      cluster = cluster
    ),
    class = c("geer_runs_test", "htest")
  )
}
