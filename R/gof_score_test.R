#' @title
#' Partition Goodness-of-Fit Tests for Binary geer Fits
#'
#' @description
#' Tests the fit of the marginal mean model of a binary response by adding
#' indicators of a partition of the observations to the linear predictor and
#' testing them with a score statistic: the deciles-of-risk test of Horton et
#' al. (1999) or the covariate-region test of Barnhart and Williamson (1998).
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param partition character string specifying the partition: groups of the
#'   fitted probabilities (\code{"risk"}) or regions of the covariates
#'   (\code{"covariates"}). See Details. Defaults to \code{"risk"}.
#' @param n_groups a single integer, at least 2, giving the number of groups of
#'   fitted probabilities for \code{partition = "risk"}. Defaults to \code{10}.
#' @param regions for \code{partition = "covariates"}, a one-sided formula
#'   whose variables, evaluated in the data used to fit the model and then in
#'   the environment of the formula, define the regions by their
#'   combinations, or a vector with one value per observation, in the order
#'   stored in \code{object}. Continuous covariates should be cut into
#'   intervals, for example with \code{\link{cut}}.
#' @param by_repeated logical indicating whether, for
#'   \code{partition = "covariates"}, the alternative model also contains
#'   indicators of the \code{repeated} positions and their interactions with
#'   the regions, as in Barnhart and Williamson (1998). Defaults to
#'   \code{TRUE}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   of the score statistic: one of the estimators described in
#'   \code{\link{vcov.geer}}, except \code{"jackknife"}. See Details. Defaults
#'   to \code{"bias-corrected"}.
#'
#' @details
#' The alternative model adds to the linear predictor the indicators of the
#' groups or regions, of which those linearly dependent on the model matrix
#' are dropped. The test of their coefficients is the score test of
#' \code{\link{anova.geer}} with \code{test = "score"}, evaluated at the fit
#' of \code{object}, so the alternative model is not fitted. The association
#' and dispersion parameters are those of \code{object}. The statistic is
#' compared with the chi-squared distribution whose degrees of freedom are the
#' number of added indicators kept.
#'
#' With \code{partition = "risk"} (Horton et al., 1999), the fitted
#' probabilities of all observations are divided at their empirical quantiles
#' of orders \eqn{0, 1/G, \ldots, 1}{0, 1/G, ..., 1}, where \eqn{G} is
#' \code{n_groups}, into intervals closed on the right. Tied fitted
#' probabilities fall into the same group, so the groups may have unequal
#' sizes and, with few distinct fitted probabilities, there may be fewer than
#' \eqn{G} groups. Implementations of Hosmer-Lemeshow type tests differ in
#' this choice and may give different results (Hosmer et al., 1997). A
#' warning is issued when fewer than \eqn{G} groups are formed. The groups
#' are labeled by their intervals of fitted probabilities.
#'
#' With \code{partition = "covariates"} (Barnhart and Williamson, 1998), the
#' regions are the combinations of the variables of \code{regions}, so
#' \code{~ treatment + gender} and \code{~ treatment:gender} define the same
#' regions. With \code{by_repeated = TRUE}, the alternative model also
#' contains the indicators of the \code{repeated} positions and of the region
#' by position combinations, which is their model (4). The cut points should
#' be chosen before looking at the data. Following Barnhart and Williamson
#' (1998), a warning is issued when a region has fewer than 10 clusters, when
#' more than 25\% of the regions have fewer than 25 clusters, or when a
#' region, or with \code{by_repeated = TRUE} a region by position cell, has no
#' events or no nonevents.
#'
#' With \code{cov_type = "naive"}, the statistic is the model-based statistic
#' \eqn{Q} of Barnhart and Williamson (1998), with the dispersion parameter of
#' \code{object}. Barnhart and Williamson fix the dispersion parameter at 1,
#' which a fit with \code{phi_fixed = TRUE} and \code{phi_value = 1}
#' reproduces. With \code{cov_type = "robust"}, the statistic is the robust
#' score statistic of Rotnitzky and Jewell (1990) used by Horton et al.
#' (1999). This differs from the robust statistic \eqn{Q_R} of Barnhart and
#' Williamson (1998), whose variance is the Schur complement of the robust
#' middle matrix of the sandwich rather than the variance of the score for the
#' added coefficients at the fit of \code{object}. In their depression example,
#' their \eqn{Q_R} of 38.93 on 8 degrees of freedom corresponds to 33.91 here.
#' The default \code{cov_type = "bias-corrected"} is that of
#' \code{\link{anova.geer}}, and \code{cov_type = "robust"} gives the
#' statistic of Horton et al. (1999).
#'
#' The papers consider the logit link. The tests are computed with the
#' derivative \eqn{D_i = \partial\mu_i/\partial\beta^\top}{D_i = d mu_i/d
#' beta'} for any link, and \code{object} can be a \code{\link{geewa}} fit of
#' a binomial family or a \code{\link{geewa_binary}} fit. The response must
#' be binary with unit prior weights.
#'
#' The tests were developed for \code{method = "gee"}. For the other
#' estimation methods, the score of the ordinary generalized estimating
#' equations is evaluated at the fitted values of \code{object}, as in
#' \code{\link{anova.geer}}, and a warning is issued, since the null
#' distribution of the statistic has not been studied for these estimators.
#'
#' @return
#' An object of class \code{c("geer_htest", "htest")}. The standard
#' components hold the score statistic (\code{statistic}) and its degrees of
#' freedom (\code{parameter}). In addition, the object contains:
#' \item{partition}{a factor giving the group or region of each observation,
#'   in the order stored in \code{object}.}
#' \item{table}{a data frame with one row per group or region and the
#'   numbers of observations and of clusters, and the observed and expected
#'   numbers of events.}
#' \item{cov_type}{the covariance matrix estimator used.}
#'
#' @references
#' Barnhart, H.X. and Williamson, J.M. (1998) Goodness-of-fit tests for GEE
#' modeling with binary responses. \emph{Biometrics}, \bold{54}, 720--729.
#'
#' Horton, N.J., Bebchuk, J.D., Jones, C.L., Lipsitz, S.R., Catalano, P.J.,
#' Zahner, G.E.P. and Fitzmaurice, G.M. (1999) Goodness-of-fit for GEE: an
#' example with mental health service utilization. \emph{Statistics in
#' Medicine}, \bold{18}, 213--222.
#'
#' Hosmer, D.W., Hosmer, T., le Cessie, S. and Lemeshow, S. (1997) A
#' comparison of goodness-of-fit tests for the logistic regression model.
#' \emph{Statistics in Medicine}, \bold{16}, 965--980.
#'
#' Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
#' parameters in semiparametric generalized linear models for cluster
#' correlated data. \emph{Biometrika}, \bold{77}, 485--497.
#'
#' @seealso \code{\link{gof_residual_test}}, \code{\link{gof_tests}},
#'   \code{\link{anova.geer}}, \code{\link{geewa}}, \code{\link{geewa_binary}}.
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
#'   orstr = "exchangeable"
#' )
#' gof_score_test(fit)
#' gof_score_test(fit, partition = "covariates", regions = ~treatment)
#'
#' @export
gof_score_test <- function(object,
                           partition = c("risk", "covariates"),
                           n_groups = 10L,
                           regions,
                           by_repeated = TRUE,
                           cov_type = "bias-corrected") {
  context <- "gof_score_test"
  data_name <- deparse1(substitute(object))
  object <- check_geer_object(object)
  partition <- match_gof_partition(partition[1L])
  cov_type <- check_gof_cov_type(
    match_cov_type(cov_type, default = "bias-corrected"),
    context
  )
  check_gof_by_repeated(by_repeated)
  risk <- identical(partition, "risk")
  if (risk && (!missing(regions) || !missing(by_repeated))) {
    stop(
      paste0(
        "'regions' and 'by_repeated' are used only with ",
        "partition = \"covariates\""
      ),
      call. = FALSE
    )
  }
  if (!risk && !missing(n_groups)) {
    stop("'n_groups' is used only with partition = \"risk\"", call. = FALSE)
  }
  if (!risk && missing(regions)) {
    stop(
      "'regions' must be supplied for partition = \"covariates\"",
      call. = FALSE
    )
  }
  n_groups <- check_integer_at_least(n_groups, "n_groups", lower = 2L)
  check_gof_fit(object, context)
  if (risk) {
    partition_factor <- build_gof_risk_groups(object, n_groups, context)
    by_repeated <- FALSE
    method <- paste(
      "Horton et al. goodness-of-fit test",
      "(groups of fitted probabilities)"
    )
  } else {
    partition_factor <- build_gof_regions(
      object,
      regions,
      env = parent.frame(),
      context = context
    )
    check_gof_regions(object, partition_factor, by_repeated)
    method <- "Barnhart-Williamson goodness-of-fit test (covariate regions)"
  }
  result <- compute_gof_score_test(
    object,
    partition_factor,
    by_repeated = by_repeated,
    cov_type = cov_type,
    context = context
  )
  structure(
    list(
      statistic = c("X-squared" = result$test_stat),
      parameter = c(df = result$test_df),
      p.value = result$test_p,
      method = sprintf("%s with the %s covariance matrix", method, cov_type),
      data.name = data_name,
      partition = partition_factor,
      table = build_gof_table(object, partition_factor),
      cov_type = cov_type
    ),
    class = c("geer_htest", "htest")
  )
}
