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
#'   original data cannot be recovered from the fitted object.
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
#' warning is issued and the asymptotic chi-squared reference is used instead;
#' this indicates a numerical problem rather than an inapplicable reference.
#' For all other missing-data patterns the large-sample chi-squared reference
#' is used silently, with no warning. Little also derives a more general
#' small-sample distribution for monotone patterns as a sum of transformed
#' independent F variables; that nonstandard reference distribution is not
#' evaluated here.
#'
#' For a \code{geer} fit, only the repeated response is tested. The function
#' reconstructs a subject-by-repeated-measure matrix from \code{id} and
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
#' Little, R.J.A. (1988) A test of missing completely at random for
#' multivariate data with missing values. \emph{Journal of the American
#' Statistical Association}, \bold{83}, 1198--1202.
#' \doi{10.1080/01621459.1988.10478722}
#'
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
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
  if (!is.numeric(maxit) || length(maxit) != 1L || is.na(maxit) ||
      !is.finite(maxit) || maxit < 1 || maxit > .Machine$integer.max ||
      maxit != floor(maxit)) {
    stop("'maxit' must be a positive integer", call. = FALSE)
  }
  maxit <- as.integer(maxit)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single positive finite number", call. = FALSE)
  }
  reference <- match.arg(reference)

  data_name <- deparse(substitute(object))
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

  x <- validate_mcar_little_matrix(x)
  if (anyNA(x) && mcar_little_has_binary_variable(x)) {
    warning(
      paste0(
        "Little (1988) notes that this MCAR test is most appropriate for ",
        "quantitative variables; consider contingency-table methods for ",
        "categorical variables"
      ),
      call. = FALSE
    )
  }

  result <- mcar_little_calculate(x, maxit = maxit, tol = tol)
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


#' @title
#' Ridout-Style Regression Diagnostic for MCAR in Longitudinal Data
#'
#' @description
#' Assesses whether longitudinal response missingness is associated with the
#' immediately preceding response, when observed, and with observed covariates. The
#' binary missingness model is fitted with \code{\link{geewa_binary}} using a
#' logit link. Occasion effects are included as nuisance terms. The same five
#' hypothesis-testing procedures available in \code{anova.geer} can be used to
#' assess response-history dependence, covariate dependence, and their joint
#' contribution.
#'
#' @param object a fitted object of class \code{"geer"}.
#' @param formula optional one-sided formula specifying fully observed
#'   covariates of response missingness, for example
#'   \code{~ treatment + age}. The previous response is added automatically
#'   and should not be included in \code{formula}. If omitted, the right-hand
#'   side of the fitted GEE mean model is used. Terms aliased with the
#'   occasion effects, such as a main linear time effect, are omitted
#'   automatically and reported in \code{dropped_covariates}.
#' @param data optional original data used to fit \code{object}. This is only
#'   needed when the original data cannot be recovered from the fitted object.
#' @param orstr working odds-ratio structure for the binary missingness GEE.
#'   One of \code{"independence"}, \code{"exchangeable"}, or
#'   \code{"unstructured"}. Defaults to \code{"independence"}.
#' @param test hypothesis-testing procedure. One of \code{"wald"},
#'   \code{"score"}, \code{"working-wald"}, \code{"working-score"}, or
#'   \code{"working-lrt"}. These are the same procedures available in
#'   \code{\link{anova.geer}}. Defaults to \code{"wald"}.
#' @param cov_type covariance estimator used for inference and for the
#'   coefficient table. One of \code{"bias-corrected"}, \code{"robust"},
#'   \code{"df-adjusted"}, \code{"jackknife"}, or \code{"naive"}. Defaults to
#'   \code{"bias-corrected"}.
#' @param pmethod approximation used to compute p-values for the modified
#'   working tests. One of \code{"rao-scott"} or
#'   \code{"satterthwaite"}. It is ignored for \code{"wald"} and
#'   \code{"score"}. Defaults to \code{"rao-scott"}.
#' @param control a list returned by \code{\link{geer_control}} controlling
#'   the binary GEE fit. Defaults to \code{geer_control()}.
#'
#' @details
#' For each subject, the diagnostic considers transitions from occasion
#' \eqn{t-1} to occasion \eqn{t}. A transition enters the risk set only when
#' \eqn{Y_{i,t-1}} is observed. The binary response is
#' \deqn{R_{it}=I(Y_{it}\text{ is missing}),}
#' and the fitted mean model has the form
#' \deqn{\mathrm{logit}\{P(R_{it}=1)\}
#' = \alpha_t + X_{it}^T\gamma + \delta Y_{i,t-1}.}
#' The occasion-specific effects \eqn{\alpha_t} are nuisance parameters.
#' The model is fitted with \code{geewa_binary(..., link = "logit",
#' method = "gee")} using the original subject identifiers.
#' If a score-based procedure is combined with \code{cov_type = "jackknife"},
#' the score and model-based information are evaluated under the corresponding
#' null model, while the covariance component is obtained from full
#' leave-one-cluster refits of the larger missingness model.
#'
#' Three tests using the procedure selected by \code{test} are returned in
#' \code{tests}:
#' \itemize{
#'   \item \code{response_history}: tests \eqn{H_0:\delta=0}. Rejection
#'   indicates that missingness depends on the immediately preceding response and
#'   is evidence against covariate-dependent MCAR/random dropout.
#'   \item \code{covariates}: tests \eqn{H_0:\gamma=0}. Covariate
#'   dependence alone is evidence against strict MCAR relative to the included
#'   covariates but can remain compatible with covariate-dependent MCAR.
#'   \item \code{overall}: tests \eqn{H_0:\gamma=0,\delta=0}.
#' }
#' The primary \code{"htest"} statistic is the \code{response_history}
#' test because this directly assesses covariate-dependent MCAR/random dropout,
#' the condition most relevant to ordinary GEE consistency. The selected
#' procedure is applied consistently to all three hypotheses. The modified
#' working LRT requires \code{orstr = "independence"}, matching the restriction
#' in \code{anova.geer}.
#'
#' The approach is closely related to Ridout's logistic-regression formulation
#' for studying random dropout. Using GEE for the binary indicators allows
#' correlation among repeated missingness indicators within a subject. The
#' procedure is diagnostic: failure to reject does not prove MCAR and does not
#' rule out dependence on unobserved responses (MNAR).
#'
#' With monotone dropout, the risk-set construction contributes each subject up
#' to the first missing response. With intermittent missingness, transitions are
#' included whenever the immediately preceding response is observed. A warning is
#' issued in that case because the result should be interpreted as a local
#' missingness-transition diagnostic rather than a pure dropout test.
#'
#' Covariates in \code{formula} must be fully observed on the reconstructed
#' longitudinal data. Completely absent measurement rows cannot be detected
#' unless they are explicitly represented in the supplied data.
#'
#' @return An object of class \code{"htest"}. In addition to the usual
#' components, it contains:
#' \itemize{
#'   \item \code{tests}: data frame containing the response-history,
#'   covariate, and overall tests using the selected procedure.
#'   \item \code{coefficients}: coefficient table for the binary GEE
#'   missingness model, or \code{NULL} when no missing response occurs in the
#'   risk set.
#'   \item \code{model}: fitted \code{geer} object returned by
#'   \code{geewa_binary}, or \code{NULL} when no missing response occurs in
#'   the risk set.
#'   \item \code{formula}: one-sided covariate formula.
#'   \item \code{dropped_covariates}: model-matrix columns omitted because
#'   they are aliased with occasion effects or earlier covariate columns.
#'   \item \code{missing}, \code{observed}, and \code{transitions}: numbers
#'   of missing outcomes, observed outcomes, and total transitions in the risk
#'   set.
#'   \item \code{clusters}: number of subjects represented in the risk set.
#'   \item \code{intermittent}: whether an observed response occurs after a
#'   missing response for at least one subject.
#'   \item \code{orstr}, \code{test}, \code{cov_type}, and
#'   \code{pmethod}: association, testing, covariance, and working-test
#'   approximation choices used for inference.
#' }
#'
#' @references
#' Rubin, D.B. (1976) Inference and missing data. \emph{Biometrika},
#' \bold{63}, 581--592. \doi{10.1093/biomet/63.3.581}
#'
#' Ridout, M.S. (1991) Testing for random dropouts in repeated measurement
#' data. \emph{Biometrics}, \bold{47}, 1617--1619.
#' \doi{10.2307/2532413}
#'
#' Fitzmaurice, G.M., Heath, A.F. and Clifford, P. (1996). Logistic regression
#' models for binary panel data with attrition. \emph{Journal of the Royal
#' Statistical Society: Series A}, \bold{159}, 249--263.
#' \doi{10.2307/2983172}
#'
#' @seealso \code{\link{mcar_little_test}}, \code{\link{geewa_binary}},
#'   \code{\link{runs_test}}.
#'
#' @examples
#' set.seed(1)
#' id <- rep(seq_len(60), each = 4)
#' time <- rep(seq_len(4), times = 60)
#' treatment <- rep(rep(c(0, 1), each = 30), each = 4)
#' y_complete <- 2 + 0.4 * treatment + 0.2 * time + rnorm(length(id))
#' y <- y_complete
#' for (i in seq_len(60)) {
#'   idx <- which(id == i)
#'   for (j in 2:4) {
#'     p <- plogis(-2.8 + 0.5 * treatment[idx[j]] +
#'       0.45 * y_complete[idx[j - 1]])
#'     if (rbinom(1, 1, p) == 1) {
#'       y[idx[j:4]] <- NA_real_
#'       break
#'     }
#'   }
#' }
#' dat <- data.frame(id, time, treatment, y)
#'
#' fit <- geewa(
#'   y ~ treatment + time,
#'   data = dat,
#'   id = id,
#'   repeated = time,
#'   family = gaussian(),
#'   corstr = "independence"
#' )
#' out <- mcar_logistic_test(fit, formula = ~ treatment, test = "wald")
#' out
#' out$tests
#'
#' ## Generalized score version of the same diagnostic
#' mcar_logistic_test(fit, formula = ~ treatment, test = "score")
#'
#' @export
mcar_logistic_test <- function(object,
                               formula = NULL,
                               data = NULL,
                               orstr = c(
                                 "independence", "exchangeable", "unstructured"
                               ),
                               test = geer_test_choices,
                               cov_type = geer_cov_type_choices,
                               pmethod = geer_pmethod_choices,
                               control = geer_control()) {
  object <- check_geer_object(object)
  orstr <- match.arg(orstr)
  opts <- normalize_geer_test_options(
    test = test[1L],
    cov_type = cov_type[1L],
    pmethod = pmethod[1L]
  )
  test <- opts$test
  cov_type <- opts$cov_type
  pmethod <- opts$pmethod
  if (identical(test, "working-lrt") && !identical(orstr, "independence")) {
    stop(
      "the modified working LRT can only be applied to an independence working model",
      call. = FALSE
    )
  }
  covariate_formula <- make_mcar_covariate_formula(object, formula)
  prepared <- reconstruct_mcar_transition_frame(
    object,
    formula = covariate_formula,
    data = data
  )

  if (prepared$intermittent) {
    warning(
      paste0(
        "intermittent response missingness detected; the diagnostic uses only ",
        "transitions for which the immediately previous response is observed ",
        "and should be interpreted as a local missingness-transition diagnostic ",
        "rather than a pure dropout test"
      ),
      call. = FALSE
    )
  }

  missing_no <- as.integer(sum(prepared$missing == 1L))
  observed_no <- as.integer(sum(prepared$missing == 0L))
  transition_no <- as.integer(length(prepared$missing))
  cluster_no <- as.integer(length(unique(prepared$id)))

  empty_tests <- data.frame(
    test = c("response_history", "covariates", "overall"),
    procedure = rep(test, 3L),
    statistic = c(0, 0, 0),
    df = c(0, 0, 0),
    p.value = c(1, 1, 1),
    row.names = NULL,
    check.names = FALSE
  )

  if (missing_no == 0L) {
    return(structure(
      list(
        statistic = c("X-squared" = 0),
        parameter = c(df = 0L),
        p.value = 1,
        method = paste0(
          "Ridout-style response-history diagnostic for MCAR using binary GEE: ",
          format_test_label(test), " test"
        ),
        data.name = "longitudinal response-missingness transitions from fitted geer object",
        alternative = "missingness depends on the previous observed response after adjustment for covariates and occasion",
        tests = empty_tests,
        coefficients = NULL,
        model = NULL,
        formula = covariate_formula,
        dropped_covariates = character(0),
        missing = missing_no,
        observed = observed_no,
        transitions = transition_no,
        clusters = cluster_no,
        intermittent = prepared$intermittent,
        orstr = orstr,
        test = test,
        cov_type = cov_type,
        pmethod = pmethod
      ),
      class = "htest"
    ))
  }
  if (observed_no == 0L) {
    stop("all risk-set outcomes are missing; the missingness model cannot be fitted", call. = FALSE)
  }
  if (cluster_no < 2L) {
    stop("the missingness model requires at least two independent clusters", call. = FALSE)
  }

  selected <- select_mcar_covariates(
    prepared$covariates,
    prepared$occasion,
    prepared$lag_response
  )

  predictor_matrix <- selected$covariates
  safe_names <- if (ncol(predictor_matrix) > 0L) {
    paste0("mcar_x", seq_len(ncol(predictor_matrix)))
  } else {
    character(0)
  }
  analysis_data <- data.frame(
    mcar_missing = prepared$missing,
    mcar_id = prepared$id,
    mcar_repeated = prepared$repeated,
    mcar_occasion = selected$occasion_factor,
    predictor_matrix,
    mcar_lag_response = prepared$lag_response,
    check.names = FALSE
  )
  if (length(safe_names)) {
    start <- 5L
    finish <- start + length(safe_names) - 1L
    names(analysis_data)[start:finish] <- safe_names
  }

  fit_terms <- c(
    if (nlevels(selected$occasion_factor) > 1L) "mcar_occasion" else character(0),
    safe_names,
    "mcar_lag_response"
  )
  fit_formula <- stats::reformulate(fit_terms, response = "mcar_missing")

  missing_fit <- fit_mcar_binary_model(
    fit_formula, analysis_data, orstr, control
  )

  nuisance_terms <- if (nlevels(selected$occasion_factor) > 1L) {
    "mcar_occasion"
  } else {
    character(0)
  }
  response_null_formula <- stats::reformulate(
    c(nuisance_terms, safe_names),
    response = "mcar_missing"
  )
  overall_null_formula <- stats::reformulate(
    nuisance_terms,
    response = "mcar_missing"
  )
  response_null_fit <- fit_mcar_binary_model(
    response_null_formula, analysis_data, orstr, control
  )
  overall_null_fit <- fit_mcar_binary_model(
    overall_null_formula, analysis_data, orstr, control
  )

  response_test <- mcar_nested_geer_test(
    response_null_fit,
    missing_fit,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )
  overall_test <- mcar_nested_geer_test(
    overall_null_fit,
    missing_fit,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )

  if (length(safe_names)) {
    covariate_null_formula <- stats::reformulate(
      c(nuisance_terms, "mcar_lag_response"),
      response = "mcar_missing"
    )
    covariate_null_fit <- fit_mcar_binary_model(
      covariate_null_formula, analysis_data, orstr, control
    )
    covariate_test <- mcar_nested_geer_test(
      covariate_null_fit,
      missing_fit,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod
    )
  } else {
    covariate_null_fit <- NULL
    covariate_test <- mcar_zero_test()
  }

  tests <- data.frame(
    test = c("response_history", "covariates", "overall"),
    procedure = rep(test, 3L),
    statistic = c(
      response_test$statistic,
      covariate_test$statistic,
      overall_test$statistic
    ),
    df = c(response_test$df, covariate_test$df, overall_test$df),
    p.value = c(
      response_test$p_value,
      covariate_test$p_value,
      overall_test$p_value
    ),
    row.names = NULL,
    check.names = FALSE
  )

  beta <- stats::coef(missing_fit)
  covariance <- stats::vcov(missing_fit, cov_type = cov_type)

  covariate_map <- if (length(safe_names)) {
    stats::setNames(selected$kept_names, safe_names)
  } else {
    character(0)
  }
  coefficient_table <- mcar_coefficient_table(
    missing_fit,
    covariance,
    covariate_map
  )

  structure(
    list(
      statistic = c("X-squared" = response_test$statistic),
      parameter = c(df = response_test$df),
      p.value = response_test$p_value,
      method = paste0(
        "Ridout-style response-history diagnostic for MCAR using binary GEE: ",
        format_test_label(test), " test"
      ),
      data.name = "longitudinal response-missingness transitions from fitted geer object",
      alternative = "missingness depends on the previous observed response after adjustment for covariates and occasion",
      tests = tests,
      coefficients = coefficient_table,
      model = missing_fit,
      formula = covariate_formula,
      dropped_covariates = selected$dropped_names,
      missing = missing_no,
      observed = observed_no,
      transitions = transition_no,
      clusters = cluster_no,
      intermittent = prepared$intermittent,
      orstr = orstr,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod
    ),
    class = "htest"
  )
}


#' @title
#' Jamshidian-Jalal Screening Diagnostic for MCAR
#'
#' @description
#' Performs the covariance-homogeneity diagnostics proposed by Jamshidian and
#' Jalal (2010) for incomplete multivariate data. The procedure groups cases by
#' their original missingness patterns, imputes the incomplete values, and then
#' assesses whether the covariance structure is homogeneous across the pattern
#' groups. It is a screening diagnostic for MCAR and does not condition on the
#' regression model fitted by \code{geer}.
#'
#' @param object A fitted \code{geer} object, a numeric matrix, or a numeric data
#'   frame containing missing values. For a fitted \code{geer} object, the
#'   original response measurements are reconstructed in subject-by-occasion
#'   form before rows omitted by the model fit are removed.
#' @param data Optional original data used to fit \code{object}. This is needed
#'   only when the data cannot be recovered from the fitted object.
#' @param method Character string selecting the diagnostic. \code{"hawkins"}
#'   uses the modified Hawkins normal-theory test, \code{"nonparametric"} uses
#'   the k-sample Anderson-Darling test, and \code{"auto"} follows the
#'   diagnostic logic of Jamshidian and Jalal: the Hawkins test is examined
#'   first and the nonparametric test distinguishes nonnormality from covariance
#'   heterogeneity when Hawkins rejects. Defaults to \code{"auto"}.
#' @param imputation Imputation method used before the diagnostics are
#'   calculated. \code{"distribution-free"} implements the residual-resampling
#'   procedure of Jamshidian and Jalal. It requires at least 10 complete cases
#'   and at least \code{2 * p} complete cases, where \code{p} is the number of
#'   variables; otherwise the function warns and uses \code{"normal"}
#'   imputation. \code{"normal"} draws from the conditional multivariate normal
#'   distribution using maximum likelihood estimates obtained by EM. Defaults
#'   to \code{"distribution-free"}.
#' @param min_pattern_size Integer greater than or equal to 2 specifying the
#'   minimum number of cases required for a missingness pattern to be retained.
#'   Defaults to \code{7}, corresponding to omitting patterns with six or fewer
#'   cases, as in the simulations and default implementation associated with
#'   Jamshidian and Jalal's procedure.
#' @param nrep Positive integer giving the number of simulated uniform samples
#'   used to obtain a p-value for a pattern-specific Neyman smooth statistic
#'   when the pattern contains fewer than \code{n_min} cases. Defaults to
#'   \code{10000}.
#' @param n_min Integer greater than or equal to 2 specifying the pattern size
#'   from which the chi-squared approximation with four degrees of freedom is
#'   used for the Neyman smooth statistic. Defaults to \code{30}.
#' @param alpha A single number strictly between 0 and 1 specifying the
#'   significance level used by the automatic diagnostic rule and its
#'   interpretation. With \code{method = "auto"}, the nonparametric diagnostic
#'   is selected when the modified Hawkins test has p-value less than or equal
#'   to \code{alpha}. Defaults to \code{0.05}.
#' @param seed \code{NULL} or a single finite number specifying the random-number
#'   seed used for imputation and simulated Neyman p-values. The previous R
#'   random-number state is restored on exit. Use \code{NULL} to use the current
#'   random-number stream. Defaults to \code{110}.
#' @param maxit Positive integer giving the maximum number of EM iterations for
#'   normal-theory imputation. Defaults to \code{1000}.
#' @param tol A single positive finite number specifying the relative
#'   convergence tolerance for the EM algorithm. Defaults to \code{1e-8}.
#'
#' @details
#' The diagnostic requires missing values. After applying
#' \code{min_pattern_size}, at least two missingness patterns must remain and at
#' least one retained pattern must contain missing values.
#'
#' The modified Hawkins component transforms within-pattern Mahalanobis
#' distances to variables that should be Uniform(0, 1) under multivariate
#' normality and homogeneous covariance matrices. A fourth-order Neyman smooth
#' test is applied within each pattern and the pattern-specific p-values are
#' combined by Fisher's method.
#'
#' The nonparametric component compares the distributions of the Hawkins
#' transformed distances across missingness-pattern groups using the k-sample
#' Anderson-Darling statistic of Scholz and Stephens (1987). This component is
#' intended to reduce sensitivity to departures from multivariate normality.
#'
#' With \code{method = "auto"}, both components are calculated. If the Hawkins
#' test does not reject, there is no evidence against its joint normality and
#' homoscedasticity null. If Hawkins rejects but the nonparametric test does not,
#' the result is consistent with nonnormality rather than covariance
#' heterogeneity. Rejection by the nonparametric test provides evidence against
#' covariance homogeneity and therefore against MCAR in the Jamshidian-Jalal
#' framework.
#'
#' The procedure assumes that, apart from the missingness-pattern grouping, the
#' cases come from a common population. Known substantive groups with genuinely
#' different covariance matrices can therefore cause rejection even when
#' missingness is MCAR. Small retained pattern groups can also make the
#' nonparametric approximation less reliable; the default
#' \code{min_pattern_size = 7} follows the practical threshold used by the
#' associated \pkg{MissMech} implementation.
#'
#' A single completed data set is used for the reported test. Jamshidian and
#' Jalal also describe multiple imputation as a sensitivity and diagnostic tool
#' for assessing imputation variability and identifying patterns that contribute
#' to a rejection. Such multiple-imputation diagnostics are not pooled by this
#' function.
#'
#' This procedure ignores the marginal regression structure and should be used
#' as a screening diagnostic. In longitudinal \code{geer} analyses,
#' \code{\link{mcar_logistic_test}} provides the complementary model-aware
#' diagnostic of whether observed covariates or previous responses predict
#' missingness. Failure to reject either component does not prove MCAR.
#'
#' @return An object inheriting from \code{"htest"} with the primary statistic
#'   and p-value for the selected diagnostic, plus the following additional
#'   components:
#'   \describe{
#'     \item{tests}{A data frame containing the Hawkins and/or nonparametric
#'       statistics and p-values.}
#'     \item{hawkins}{Detailed results from the modified Hawkins test.}
#'     \item{nonparametric}{Detailed results from the k-sample
#'       Anderson-Darling test, when calculated.}
#'     \item{patterns}{The retained missingness patterns, with 1 denoting a
#'       missing value and 0 an observed value.}
#'     \item{pattern.counts}{Numbers of cases in the retained patterns.}
#'     \item{omitted.patterns}{Patterns omitted because they did not meet
#'       \code{min_pattern_size}.}
#'     \item{n}{Number of rows retained for the diagnostic.}
#'     \item{p}{Number of variables or repeated measurements tested.}
#'     \item{imputation.requested,imputation.used}{Requested and actual
#'       imputation methods.}
#'     \item{complete.cases}{Number of complete cases used by the imputation
#'       procedure.}
#'     \item{imputed.data}{The single completed data set used by the tests.}
#'     \item{location}{Estimated location vector used for imputation.}
#'     \item{covariance}{Estimated covariance matrix used for imputation.}
#'     \item{em.iterations}{Number of EM iterations used for normal-theory
#'       imputation.}
#'     \item{em.converged}{Whether the EM algorithm converged for normal-theory
#'       imputation.}
#'     \item{method.requested}{Diagnostic requested through \code{method}.}
#'     \item{selected.test}{Diagnostic supplying the primary statistic and
#'       p-value.}
#'     \item{alpha}{Significance level used by the automatic diagnostic rule.}
#'     \item{interpretation}{A concise interpretation using \code{alpha}.}
#'   }
#'
#' @references
#' Jamshidian, M. and Jalal, S. (2010) Tests of homoscedasticity, normality,
#' and missing completely at random for incomplete multivariate data.
#' \emph{Psychometrika}, \bold{75}, 649--674. \doi{10.1007/s11336-010-9175-3}.
#'
#' Jamshidian, M., Jalal, S. and Jansen, C. (2014) MissMech: An R package for
#' testing homoscedasticity, multivariate normality, and missing completely at
#' random (MCAR). \emph{Journal of Statistical Software}, \bold{56}, 1--31.
#' \doi{10.18637/jss.v056.i06}.
#'
#' Hawkins, D.M. (1981) A new test for multivariate normality and
#' homoscedasticity. \emph{Technometrics}, 23, 105--110.
#'
#' Scholz, F. W. and Stephens, M. A. (1987). K-sample Anderson-Darling tests.
#' \emph{Journal of the American Statistical Association}, 82, 918--924.
#'
#' @seealso \code{\link{mcar_little_test}}, \code{\link{mcar_logistic_test}}
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(240), ncol = 3)
#' x[41:60, 3] <- NA
#' x[61:80, 2:3] <- NA
#'
#' mcar_homoscedasticity_test(
#'   x,
#'   method = "nonparametric",
#'   nrep = 100
#' )
#'
#' @export
mcar_homoscedasticity_test <- function(
    object,
    data = NULL,
    method = c("auto", "nonparametric", "hawkins"),
    imputation = c("distribution-free", "normal"),
    min_pattern_size = 7L,
    nrep = 10000L,
    n_min = 30L,
    alpha = 0.05,
    seed = 110L,
    maxit = 1000L,
    tol = 1e-8) {
  method <- match.arg(method)
  imputation <- match.arg(imputation)

  if (length(min_pattern_size) != 1L || is.na(min_pattern_size) ||
      min_pattern_size < 2 || min_pattern_size != as.integer(min_pattern_size)) {
    stop("'min_pattern_size' must be an integer greater than or equal to 2", call. = FALSE)
  }
  min_pattern_size <- as.integer(min_pattern_size)

  if (length(nrep) != 1L || is.na(nrep) || nrep < 1 || nrep != as.integer(nrep)) {
    stop("'nrep' must be a positive integer", call. = FALSE)
  }
  nrep <- as.integer(nrep)

  if (length(n_min) != 1L || is.na(n_min) || n_min < 2 || n_min != as.integer(n_min)) {
    stop("'n_min' must be an integer greater than or equal to 2", call. = FALSE)
  }
  n_min <- as.integer(n_min)

  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number strictly between 0 and 1", call. = FALSE)
  }
  if (length(maxit) != 1L || is.na(maxit) || maxit < 1 || maxit != as.integer(maxit)) {
    stop("'maxit' must be a positive integer", call. = FALSE)
  }
  maxit <- as.integer(maxit)
  if (length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single positive finite number", call. = FALSE)
  }

  is_geer <- inherits(object, "geer")
  x <- if (is_geer) {
    extract_geer_mcar_matrix(object, data = data)
  } else {
    object
  }
  x <- validate_mcar_homoscedasticity_matrix(x)
  if (!anyNA(x)) {
    stop(
      "the Jamshidian-Jalal diagnostic requires missing values",
      call. = FALSE
    )
  }

  pattern_info <- jj_pattern_information(x, min_pattern_size = min_pattern_size)
  x_used <- validate_mcar_homoscedasticity_matrix(pattern_info$x)
  if (!anyNA(x_used)) {
    stop(
      "the Jamshidian-Jalal diagnostic requires missing values in at least one retained pattern",
      call. = FALSE
    )
  }

  data_name <- deparse1(substitute(object))

  result <- jj_with_seed(seed, {
    imputed <- jj_impute(
      x_used,
      imputation = imputation,
      maxit = maxit,
      tol = tol
    )
    if (anyNA(imputed$data) || any(!is.finite(imputed$data))) {
      stop("imputation produced missing or non-finite completed values", call. = FALSE)
    }

    hawkins <- jj_hawkins_test(
      completed = imputed$data,
      group = pattern_info$group,
      group_counts = pattern_info$group_counts,
      nrep = nrep,
      n_min = n_min,
      test_uniformity = !identical(method, "nonparametric")
    )

    nonparametric <- NULL
    if (method %in% c("auto", "nonparametric")) {
      nonparametric <- jj_nonparametric_test(
        hawkins = hawkins,
        group = pattern_info$group,
        group_counts = pattern_info$group_counts
      )
    }

    list(imputed = imputed, hawkins = hawkins, nonparametric = nonparametric)
  })

  if (identical(method, "hawkins")) {
    statistic <- c("Fisher chi-squared" = result$hawkins$statistic)
    parameter <- c(df = result$hawkins$parameter)
    p_value <- result$hawkins$p.value
    method_label <- "Jamshidian-Jalal modified Hawkins MCAR screening diagnostic"
    selected_test <- "hawkins"
  } else if (identical(method, "nonparametric")) {
    statistic <- c("Anderson-Darling" = result$nonparametric$statistic)
    parameter <- NULL
    p_value <- result$nonparametric$p.value
    method_label <- "Jamshidian-Jalal nonparametric MCAR screening diagnostic"
    selected_test <- "nonparametric"
  } else if (result$hawkins$p.value > alpha) {
    statistic <- c("Fisher chi-squared" = result$hawkins$statistic)
    parameter <- c(df = result$hawkins$parameter)
    p_value <- result$hawkins$p.value
    method_label <- "Jamshidian-Jalal automatic MCAR screening diagnostic (modified Hawkins)"
    selected_test <- "hawkins"
  } else {
    statistic <- c("Anderson-Darling" = result$nonparametric$statistic)
    parameter <- NULL
    p_value <- result$nonparametric$p.value
    method_label <- "Jamshidian-Jalal automatic MCAR screening diagnostic (nonparametric)"
    selected_test <- "nonparametric"
  }

  tests <- data.frame(
    test = character(),
    statistic = numeric(),
    df = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )
  if (!identical(method, "nonparametric")) {
    tests <- rbind(
      tests,
      data.frame(
        test = "hawkins",
        statistic = result$hawkins$statistic,
        df = result$hawkins$parameter,
        p.value = result$hawkins$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
  if (!is.null(result$nonparametric)) {
    tests <- rbind(
      tests,
      data.frame(
        test = "nonparametric",
        statistic = result$nonparametric$statistic,
        df = NA_real_,
        p.value = result$nonparametric$p.value,
        stringsAsFactors = FALSE
      )
    )
  }

  out <- list(
    statistic = statistic,
    parameter = parameter,
    p.value = p_value,
    method = method_label,
    data.name = data_name,
    alternative = if (identical(selected_test, "hawkins")) {
      "multivariate normality or covariance homogeneity fails"
    } else {
      "covariance structure differs across missingness-pattern groups"
    },
    tests = tests,
    hawkins = result$hawkins,
    nonparametric = result$nonparametric,
    patterns = pattern_info$pattern_matrix,
    pattern.counts = stats::setNames(
      pattern_info$group_counts,
      rownames(pattern_info$pattern_matrix)
    ),
    omitted.patterns = pattern_info$omitted_patterns,
    n = nrow(x_used),
    p = ncol(x_used),
    imputation.requested = result$imputed$requested,
    imputation.used = result$imputed$used,
    complete.cases = result$imputed$n_complete,
    imputed.data = result$imputed$data,
    location = result$imputed$mu,
    covariance = result$imputed$sigma,
    em.iterations = result$imputed$iterations,
    em.converged = result$imputed$converged,
    method.requested = method,
    selected.test = selected_test,
    alpha = alpha,
    interpretation = jj_interpret(
      method = method,
      hawkins = result$hawkins,
      nonparametric = result$nonparametric,
      alpha = alpha
    )
  )
  class(out) <- c("mcar_homoscedasticity_test", "htest")
  out
}
