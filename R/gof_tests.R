#' @title
#' Goodness-of-Fit Tests for Binary geer Fits
#'
#' @description
#' Computes the goodness-of-fit tests of \code{\link{gof_score_test}} and
#' \code{\link{gof_residual_test}} for a binary \code{geer} fit and collects
#' them in one table.
#'
#' @inheritParams gof_score_test
#' @inheritParams gof_residual_test
#' @param n_groups a single integer, at least 2, giving the number of groups of
#'   fitted probabilities for the test of Horton et al. (1999). Defaults to
#'   \code{10}.
#' @param regions optional one-sided formula or vector defining covariate
#'   regions, as in \code{\link{gof_score_test}}. If supplied, the
#'   covariate-region tests of Barnhart and Williamson (1998) are included.
#' @param by_repeated logical indicating whether the alternative model of the
#'   covariate-region tests also contains indicators of the \code{repeated}
#'   positions and their interactions with the regions, as in
#'   \code{\link{gof_score_test}}. Used only when \code{regions} is supplied.
#'   Defaults to \code{TRUE}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   of the score statistics, as in \code{\link{gof_score_test}}. Defaults to
#'   \code{"bias-corrected"}.
#'
#' @details
#' The table contains the deciles-of-risk test of Horton et al. (1999) with
#' \code{cov_type}, the covariate-region tests of Barnhart and Williamson
#' (1998) with the naive covariance matrix and with \code{cov_type} when
#' \code{regions} is supplied, and the Pearson and unweighted sum of squares
#' tests of Pan (2002) with \code{covariance}. The tests have power against
#' different departures from the model (Evans and Li, 2005), so they
#' complement each other, and a single non-significant test is weak evidence
#' of fit.
#'
#' @return
#' An object of class \code{c("geer_anova", "anova", "data.frame")} with one
#' row per test and columns \code{Statistic}, \code{Df} (for the score
#' tests), \code{Mean} and \code{Variance} (the approximate mean and variance
#' of the residual statistics), and \code{Pr(Test)}, the p-value of the
#' chi-squared test or of the two-sided normal test.
#'
#' @references
#' Barnhart, H.X. and Williamson, J.M. (1998) Goodness-of-fit tests for GEE
#' modeling with binary responses. \emph{Biometrics}, \bold{54}, 720--729.
#'
#' Evans, S. and Li, L. (2005) A comparison of goodness of fit tests for the
#' logistic GEE model. \emph{Statistics in Medicine}, \bold{24}, 1245--1261.
#'
#' Horton, N.J., Bebchuk, J.D., Jones, C.L., Lipsitz, S.R., Catalano, P.J.,
#' Zahner, G.E.P. and Fitzmaurice, G.M. (1999) Goodness-of-fit for GEE: an
#' example with mental health service utilization. \emph{Statistics in
#' Medicine}, \bold{18}, 213--222.
#'
#' Pan, W. (2002) Goodness-of-fit tests for GEE with correlated binary data.
#' \emph{Scandinavian Journal of Statistics}, \bold{29}, 101--110.
#'
#' @seealso \code{\link{gof_score_test}}, \code{\link{gof_residual_test}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
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
#' gof_tests(fit)
#' gof_tests(fit, regions = ~treatment)
#'
#' @export
gof_tests <- function(object,
                      n_groups = 10L,
                      regions,
                      by_repeated = TRUE,
                      cov_type = "bias-corrected",
                      covariance = c("unstructured", "empirical")) {
  context <- "gof_tests"
  object <- check_geer_object(object)
  n_groups <- check_integer_at_least(n_groups, "n_groups", lower = 2L)
  cov_type <- check_gof_cov_type(
    match_cov_type(cov_type, default = "bias-corrected"),
    context
  )
  covariance <- match_gof_covariance(covariance[1L])
  check_gof_by_repeated(by_repeated)
  if (missing(regions) && !missing(by_repeated)) {
    stop(
      "'by_repeated' is used only when 'regions' is supplied",
      call. = FALSE
    )
  }
  check_gof_fit(object, context)
  score_row <- function(result) {
    c(result$test_stat, result$test_df, NA_real_, NA_real_, result$test_p)
  }
  rows <- list()
  risk <- build_gof_risk_groups(object, n_groups, context)
  rows[[sprintf("Horton et al., %s", cov_type)]] <- score_row(
    compute_gof_score_test(object, risk, FALSE, cov_type, context)
  )
  if (!missing(regions)) {
    region <- build_gof_regions(
      object,
      regions,
      env = parent.frame(),
      context = context
    )
    check_gof_regions(object, region, by_repeated)
    for (type in unique(c("naive", cov_type))) {
      rows[[sprintf("Barnhart-Williamson, %s", type)]] <- score_row(
        compute_gof_score_test(object, region, by_repeated, type, context)
      )
    }
  }
  components <- compute_gof_residual_parts(object, covariance)
  for (statistic in geer_gof_statistics) {
    result <- compute_gof_residual_test(components, statistic, context)
    label <- if (identical(statistic, "pearson")) "Pearson" else "USS"
    rows[[sprintf("Pan %s, %s", label, covariance)]] <- c(
      result$value, NA_real_, result$mean, result$variance, result$p_value
    )
  }
  ans <- as.data.frame(do.call(rbind, rows))
  names(ans) <- c("Statistic", "Df", "Mean", "Variance", "Pr(Test)")
  attr(ans, "heading") <- c(
    "Goodness-of-Fit Tests for Binary Responses\n",
    paste0(
      "Score tests: chi-squared on Df. Pan tests: two-sided normal on ",
      "(Statistic - Mean)/sqrt(Variance).\n"
    )
  )
  class(ans) <- c("geer_anova", "anova", "data.frame")
  ans
}
